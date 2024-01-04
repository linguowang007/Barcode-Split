import gzip
import os
import argparse
import logging
import time
import itertools
from collections import Counter, defaultdict
from multiprocessing import Pool
from array import array
from .dependence import samtools, tabix, mawk, bgzip
from .exception import SplitBAMError
from .split_bam_file import de_duplicates, estimate_break_position


def chunk_sort(sam_file, chunk_file, header_file, s_pos, q_len, tag='CB', tag_type='Z'):
    pre_len = len(f'{tag}:{tag_type}:') + 1  # CB:Z:
    cmd = rf'''{bgzip} -b {s_pos} -s {q_len} {sam_file} | \
            cat {header_file} - | \
            {samtools} sort -t CB -O SAM -T {chunk_file}- - | \
            {samtools} view | {bgzip} -i -I {chunk_file}.gzi > {chunk_file}
        zcat {chunk_file} | \
            {mawk} '{{
                for(i=12; i<=NF; i++){{
                    if($i~/^{tag}:{tag_type}:/){{cb=substr($i,{pre_len}); print cb"\t"++see[cb]"\t"length+1; next
                        }}
                    }}
                }}' | {bgzip} > {chunk_file}.tag.gz
        {tabix} -s1 -b2 -e2 -C {chunk_file}.tag.gz'''
    with open(f'{chunk_file}.sh', 'wt') as h1:
        print(cmd, file=h1)
    os.system(f'sh {chunk_file}.sh && rm {chunk_file}.sh')
    
    inner_counter = Counter()
    with gzip.open(f'{chunk_file}.tag.gz', 'rt') as h1:
        for recode in h1:
            parts = recode.split()
            inner_counter[parts[0]] += int(parts[-1])
    with os.popen(f'{tabix} -l {chunk_file}.tag.gz') as h1:
        inner_cell_order = [i.strip() for i in h1]
    os.system(f'rm {chunk_file}.tag.gz {chunk_file}.tag.gz.csi')
    return chunk_file, inner_counter, inner_cell_order


def second_sort(out_file, header_file, group_name, group_row, tag='CB'):
    inter_cmd = [f'cp {header_file} {out_file}']
    inter_cmd.extend([f'{bgzip} -b {s_pos} -s {q_len} {chunk_file} | bgzip >> {out_file}'
                      for chunk_file, s_pos, q_len, _ in group_row])
    inter_cmd.append(rf'''{samtools} sort -t {tag} -O SAM \
        -T {out_file}.body.gz- {out_file} | \
        {samtools} view | {bgzip} -i -I {out_file}.body.gz.gzi > {out_file}.body.gz''')
    with open(f'{out_file}.sh', 'wt') as h1:
        print('\n'.join(inter_cmd), file=h1)
    os.system(f'sh {out_file}.sh && rm {out_file}.sh {out_file}')
    return group_name, f'{out_file}.body.gz'


def split_bam_by_tag_fast(*, bam, tag_list, out_dir, nt=16, tag='CB', tag_type='Z'):
    """
            Help Document for split_bam_by_tag_fast

            Function Name: split_bam_by_tag_fast

            Description:
            The split_bam_by_tag_fast function is designed to split BAM files by specific tags.
            It's particularly fast for multiprocessing computing in large-scale single-cell sequencing data
            where differentiating data based on certain cell barcodes or tags is required.

            Parameters:
            - bam (string): Path to the input BAM file.
            - tag_list (string): Path to the file containing a list of tags.
            - out_dir (string): Path to the output directory where the split BAM files will be saved.
            - nt (int, optional): Number of threads to use for processing. Defaults to 16.
            - tag (string, optional): The specific tag to be used for splitting the BAM file. Defaults to 'CB'.
            - tag_type (string, optional): Type of the tag as defined in the BAM file. Defaults to 'Z'.

            Returns:
            A dictionary where keys are CB IDs and values are paths to the split sub-BAM files.

            Example Usage:
            split_bam_by_tag_fast(bam="path/to/bamfile.bam",
                             tag_list="path/to/tags.txt",
                             out_dir="path/to/output",
                             nt=8,
                             tag="CB",
                             tag_type="Z")
            """
    logging.basicConfig(filename=f'{out_dir}.log', filemode='w', level=logging.INFO, format='%(message)s')
    
    if os.path.exists(out_dir):
        raise SplitBAMError(f'{out_dir} already exists.')
    
    start_time = time.monotonic()
    t1 = time.monotonic()
    sam = f'{out_dir}/to-bgzip.bam.gz'
    cmd1 = rf'''mkdir {out_dir}
        {samtools} view -@ {nt} --tag-file {tag}:{tag_list} {bam} | {bgzip} -@ {nt} -i -I {sam}.gzi > {sam}'''
    os.system(cmd1)
    logging.info(f'extract and bgzip: {time.monotonic() - t1} sec.')
    
    t1 = time.monotonic()
    with open(f'{sam}.gzi', 'rb') as f, open(sam, 'rb') as f1:
        index = array('q')
        index.fromfile(f, 1)
        index.fromfile(f, index[0] * 2)
        last_index = index[2::2][-1]
    with os.popen(f'{bgzip} -b {last_index} {sam} | wc -c') as f:
        sam_size = last_index + int([i.strip() for i in f][0])
    chunk_size = int(sam_size / nt / 2) + 1
    chunk_index = estimate_break_position(bgzip_file=sam, chunk_size=chunk_size)
    logging.info(f'break positions: {time.monotonic() - t1} sec.')
    
    t1 = time.monotonic()
    header = f'{sam}.header'
    os.system(f'{samtools} view -H {bam} > {header}')
    # def chunk_sort(sam_file, chunk_file, header_file, s_pos, q_len, tag='CB', tag_type='Z'):
    args = [(sam, f'{sam}-{num}.gz', header, pos1, len1, tag, tag_type)
            for num, (pos1, len1), in enumerate(chunk_index)]
    
    with Pool(nt) as p:
        result = p.starmap(chunk_sort, args, chunksize=1)
    logging.info(f'1st sorted: {time.monotonic() - t1} sec.')
    
    t1 = time.monotonic()
    d1 = {a: (b, c) for a, b, c in result if c}  # chunk-file, counter, cell-order, cell-order not empty
    chunk_file_order = [i[0] for i in result if i[0] in d1]  # order chunk file, not empty
    final_cell_order = sorted(set(itertools.chain.from_iterable(i[-1] for i in result)))
    
    total_index = Counter()
    for c2, _ in d1.values():
        total_index.update(c2)
    
    chunk_size = int(sum(total_index.values()) / nt / 2) + 1
    acc_cel_len = 0
    group = 1
    cell_group = {}
    for cell in final_cell_order:
        if acc_cel_len > chunk_size:
            group += 1
            acc_cel_len = 0
        cell_group[cell] = f'group{group}'
        acc_cel_len += total_index[cell]
    
    final_index = []  # chunk-file, cell-barcode, start-pos, query-len, group-ID
    group_index = []  # chunk-file, start-pos, query-len, group-ID
    for chunk in chunk_file_order:
        file_count, barcode_order = d1[chunk]
        sub_cell_len = [file_count[i] for i in barcode_order]
        sub_start_pos = [0] + list(itertools.accumulate(sub_cell_len))
        sub_group_name = [cell_group[i] for i in barcode_order]
        final_index.extend([(chunk, cell_barcode, pos1, len1, sub_g)
                            for cell_barcode, pos1, len1, sub_g in
                            zip(barcode_order, sub_start_pos, sub_cell_len, sub_group_name)])
        
        order_uniq_group_id = de_duplicates(sub_group_name)
        group_counter = Counter()
        for group_id, len1 in zip(sub_group_name, sub_cell_len):
            group_counter[group_id] += len1
        group_len = [group_counter[i] for i in order_uniq_group_id]
        group_start_pos = [0] + list(itertools.accumulate(group_len))
        group_index.extend([(chunk, pos1, len1, sub_g)
                            for pos1, len1, sub_g in zip(group_start_pos, group_len, order_uniq_group_id)])
    
    logging.info(f'make group index: {time.monotonic() - t1} sec.')
    
    t1 = time.monotonic()
    final_group_index = defaultdict(list)
    for items in group_index:
        final_group_index[items[-1]].append(items)
    
    header2 = f'{sam}.header.gz'
    os.system(f'{samtools} view -H {bam} | {bgzip} > {header2} ')
    # def second_sort(out_file, header_file, group_name, group_row, tag='CB'):
    args = [(f'{sam}-{g_name}.gz', header2, g_name, g_rows) for g_name, g_rows in final_group_index.items()]
    with Pool(nt) as p:
        group_sort_result = {a: b for a, b in p.starmap(second_sort, args, chunksize=1)}  # group-ID, bgzip file
    logging.info(f'2nd sorted: {time.monotonic() - t1} sec.')
    
    t1 = time.monotonic()
    tmp = defaultdict(list)
    for items in final_index:
        tmp[items[-1]].append(items)
    
    sort_tag_index = []  # chunk file, barcode-name, start-pos, query-len
    for group_id, sorted_chunk in group_sort_result.items():
        keep_row = tmp[group_id]
        c1 = Counter()
        for sub_row in keep_row:
            c1[sub_row[1]] += sub_row[-2]  # cell-name as key, accumulate the text length
        sub_cell_order = sorted(c1)
        sub_cell_len = [c1[i] for i in sub_cell_order]
        sub_start_pos = [0] + list(itertools.accumulate(sub_cell_len))
        sort_tag_index.extend([sorted_chunk, cell_name, pos1, len1]
                              for cell_name, pos1, len1 in zip(sub_cell_order, sub_start_pos, sub_cell_len))
    
    cmds = [rf'''{bgzip} -b {pos1} -s {len1} {chunk_file} | cat {header} - | \
        {samtools} view -h --write-index -o {out_dir}/{cell_name}.bam -'''
            for chunk_file, cell_name, pos1, len1 in sort_tag_index]
    with Pool(nt) as p:
        p.map(os.system, cmds, chunksize=1)
    
    tmp1 = chunk_file_order + [f'{i}.gzi' for i in chunk_file_order]
    os.system(rf'''rm {' '.join(tmp1)}''')
    tmp2 = [i[0] for i in de_duplicates(sort_tag_index, key=lambda x: x[0])]
    tmp2 = tmp2 + [f'{i}.gzi' for i in tmp2]
    os.system(rf'''rm {' '.join(tmp2)} {header} {header2} {sam} {sam}.gzi''')
    
    logging.info(f'split tag to sub BAMs and remove tmp files: {time.monotonic() - t1} sec.')
    logging.info(f'total time:{time.monotonic() - start_time} sec.')
    return {cell_name: f'{out_dir}/{cell_name}.bam' for _, cell_name, *_ in sort_tag_index}


# Define a function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Split BAM files by Tag Name using Split-Barcode.')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--tag_list', required=True, help='List of tags')
    parser.add_argument('--out_dir', required=True, help='Output directory')
    parser.add_argument('--nt', type=int, default=16, help='Number of threads (default: 16)')
    parser.add_argument('--tag', default='CB', help="Tag name (default: 'CB')")
    parser.add_argument('--tag_type', default='Z', help="Tag values type (default: 'Z')")
    return parser.parse_args()


def main():
    args = parse_arguments()
    # bam, tag_list, out_dir, nt=16, tag='CB', tag_type
    split_bam_by_tag_fast(bam=args.bam, tag_list=args.tag_list, out_dir=args.out_dir, nt=args.nt, tag=args.tag,
                          tag_type=args.tag_type)


if __name__ == '__main__':
    main()
