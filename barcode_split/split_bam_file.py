import gzip
import os
import time
import logging
import itertools
import argparse
from collections import Counter
from .exception import SplitBAMError
from .dependence import samtools, tabix, bgzip, mawk
from multiprocessing import Pool


def pairwise(iterable):  # ABCD -> (A, B), (B, C), (C, D)
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def de_duplicates(input_list, key=None):
    # de_duplicates('1 2 4 2 5 4 6'.split()) --> 1 2 4 5 6
    seen = set()
    dedup_list = []
    for item in input_list:
        key_value = item if key is None else key(item)
        if key_value not in seen:
            seen.add(key_value)
            dedup_list.append(item)
    return dedup_list


def estimate_break_position(bgzip_file, chunk_size=1000000000):
    # chunk the bgzip file in chunks_size and ensure all chunk end with '\n'
    break_pos = 0
    break_pos_list = [break_pos]
    while True:
        break_pos += chunk_size
        with os.popen(f'{bgzip} -b {break_pos} {bgzip_file} | head -1') as f:
            try:
                extend = [len(i) for i in f][0]  # when out of file end, meet IndexError
                break_pos += extend
                break_pos_list.append(break_pos)
            except IndexError:
                break
    break_pos_list.append(break_pos)  # add the last position
    return [(a1, a2 - a1) for a1, a2 in pairwise(break_pos_list)]  # start, size


def get_chunk_index(bgzip_sam, pos1, len1, out_file, tag='CB', tag_type='Z'):
    pre_len = len(f'{tag}:{tag_type}:') + 1
    index_cmd = rf'''{bgzip} -b {pos1} -s {len1} {bgzip_sam} | \
        {mawk} '{{for (i=12; i<=NF; i++) {{
            if($i~/^{tag}:{tag_type}:/){{cb=substr($i, {pre_len}); print cb"\t"++see[cb]"\t"length+1; next;}} }}
        }}' | {bgzip}  > {out_file}
        {tabix} -s1 -b2 -e2 -C {out_file}'''
    with open(f'{out_file}.sh', 'wt') as h1:
        print(index_cmd, file=h1)
    os.system(f'sh {out_file}.sh && rm {out_file}.sh')
    in_counter = Counter()
    with gzip.open(f'{out_file}', 'rt') as h1:
        for line in h1:
            parts = line.split()
            in_counter[parts[0]] += int(parts[-1])
    with os.popen(f'{tabix} -l {out_file}') as h1:
        in_cell = [i.strip() for i in h1]
    os.system(f'rm {out_file} {out_file}.csi')
    return in_counter, in_cell


def split_bam_by_tag(*, bam, tag_list, out_dir, nt=16, tag='CB', tag_type='Z'):
    """
        Help Document for split_bam_by_tag

        Function Name: split_bam_by_tag

        Description:
        The split_bam_by_tag function is designed to split BAM files by specific tags.
        It's particularly useful in processing large-scale single-cell sequencing data
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
        split_bam_by_tag(bam="path/to/bamfile.bam",
                         tag_list="path/to/tags.txt",
                         out_dir="path/to/output",
                         nt=8,
                         tag="CB",
                         tag_type="Z")
        """
    if os.path.exists(out_dir):
        raise SplitBAMError(f'{out_dir} already exists.')
    
    logging.basicConfig(format='%(message)s', level=logging.INFO, filename=f'{out_dir}.log')
    sam_file = f'{out_dir}/sub-set.sam.gz'
    with Pool(nt) as p:
        start_time = time.monotonic()
        t1 = time.monotonic()
        # sort bam by TAG and bgzip and index
        sort_cmd = rf'''mkdir -p {out_dir}
                {samtools} view -h -@ {nt} --tag-file {tag}:{tag_list} {bam} | \
                {samtools} sort -@ {nt} -t {tag} -O SAM -T {sam_file}- - | \
                {samtools} view -@ {nt} | \
                {bgzip} -@ {nt} -i -I {sam_file}.gzi > {sam_file}'''
        os.system(sort_cmd)
        logging.info(f'Sort by {tag}: {round(time.monotonic() - t1, 2)} sec.')
        
        t1 = time.monotonic()
        break_pos_list = estimate_break_position(sam_file)  # start pos in sam file, query length
        # def get_chunk_index(sam_file, pos1, len1, out_file, tag='CB', tag_type='Z'):
        args = [(sam_file, pos1, len1, f'{sam_file}-{num}.gz', tag, tag_type) for num, (pos1, len1) in
                enumerate(break_pos_list)]
        result1 = p.starmap(get_chunk_index, args, chunksize=1)
        final_index = Counter()
        for c, _ in result1:  # counter, cell order
            final_index.update(c)
        cell_order = de_duplicates(itertools.chain.from_iterable(i[-1] for i in result1))
        cell_len = [final_index[i] for i in cell_order]
        start_pos = [0] + list(itertools.accumulate(final_index[i] for i in cell_order))
        final_index = list(zip(cell_order, start_pos, cell_len))  # barcode, start, length
        logging.info(f'Make index: {round(time.monotonic() - t1, 2)} sec.')
        
        t1 = time.monotonic()
        # split sorted BAM file by tag IDs
        header = f'{sam_file}.header'
        os.system(f'{samtools} view -H {bam} > {header}')
        cmds = [rf'''{bgzip} -b {t1} -s {t2} {sam_file} | cat {header} - | \
                    {samtools} view --write-index -o {out_dir}/{cell}.sort.bam''' for cell, t1, t2 in final_index]
        p.map(os.system, cmds, chunksize=10)
        os.system(f'rm {header} {sam_file} {sam_file}.gzi')
        
        logging.info(f'Split file: {round(time.monotonic() - t1, 2)} sec.')
    logging.info(f'Split {len(cell_order)} barcodes ({out_dir}) total: {round(time.monotonic() - start_time, 2)}s.')
    
    return {cell: f'{out_dir}/{cell}.sort.bam' for cell, *_ in final_index}


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
    split_bam_by_tag(bam=args.bam, tag_list=args.tag_list, out_dir=args.out_dir, nt=args.nt, tag=args.tag,
                     tag_type=args.tag_type)


if __name__ == '__main__':
    main()
