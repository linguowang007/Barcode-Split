import os
import time
import logging
import itertools
import argparse
from collections import defaultdict
from .exception import SplitBAMError
from .dependence import samtools, tabix, bgzip, mawk
from multiprocessing import Pool, set_start_method


def pairwise(iterable):  # ABCD -> (A, B), (B, C), (C, D)
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def get_tag_chunk(chunk_file):
    barcode_size = []
    with os.popen(f'{tabix} -l {chunk_file}') as h1:
        tokens = [i.strip() for i in h1]
    for i in tokens:
        with os.popen(f'{tabix} {chunk_file} {i}') as h1:
            barcode_size.append([i, sum(int(line.split()[-1]) for line in h1)])
    return barcode_size


def dedup(items):  # drop the duplicated value in iterable, while keeping order
    see = set()
    for i in items:
        if i not in see:
            yield i
            see.add(i)


def estimate_break_position(bgzip_file, chunk_size=1000000000):
    # chunk the bgzip file in chunks_size and ensure all chunk end with '\n'
    break_pos = 0
    break_pos_list = [break_pos]
    while True:
        break_pos += chunk_size
        with os.popen(f'bgzip -b {break_pos} {bgzip_file} | head -1 ') as f:
            try:
                extend = [len(i) for i in f][0]  # when out of file end, meet IndexError
                break_pos += extend
                break_pos_list.append(break_pos)
            except IndexError:
                break
    break_pos_list.append(break_pos)  # add the last position
    return [(a1, a2 - a1) for a1, a2 in pairwise(break_pos_list)]  # start, size


def split_bam_by_tag(bam, tag_list, out_dir, nt=16, tag='CB', tag_type='Z'):
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
    with open(tag_list) as f:
        tag_values = [i.strip() for i in f]
    if 'header' in tag_values:  # "header" was reserved word
        raise SplitBAMError('Barcode list contain "header" as values.')
    logging.basicConfig(format='%(message)s', level=logging.INFO, filename=f'{out_dir}.log')

    with Pool(nt) as p:
        start_time = time.monotonic()

        t1 = time.monotonic()
        # sort bam by TAG and bgzip and index
        sort_cmd = rf'''mkdir -p {out_dir}
                {samtools} view -h -@ {nt} --tag-file {tag}:{tag_list} {bam} | \
                {samtools} sort -@ {nt} -t {tag} -O SAM -T {out_dir}/sub-set - | \
                {bgzip} -@ {nt} -i -I {out_dir}/sub-set.sam.gz.gzi > {out_dir}/sub-set.sam.gz'''
        os.system(sort_cmd)
        if not os.path.exists(f'{out_dir}/sub-set.sam.gz.gzi'):
            raise SplitBAMError('Fail to index the sorted SAM file.')
        logging.info(f'Sort by {tag}: {round(time.monotonic() - t1, 2)} sec.')

        t1 = time.monotonic()
        # index each chunk with multiple threads
        break_pos_list = estimate_break_position(f'{out_dir}/sub-set.sam.gz')
        pre_len = len(f'{tag}:{tag_type}:') + 1
        index_cmd = [rf'''{bgzip} -b {start} -s {size} {out_dir}/sub-set.sam.gz | \
            {mawk} '/^@/{{print "header\t"NR"\t"length+1; next}}
            {{for (i=12; i<=NF; i++) {{
                if($i~/^{tag}:{tag_type}:/){{cb=substr($i, {pre_len}); print cb"\t"see[cb]++"\t"length+1; next;}} }}
            }}' | {bgzip}  > {out_dir}/sub-set-{num}.sam.tag.gz
            {tabix} -s1 -b2 -e2 -C {out_dir}/sub-set-{num}.sam.tag.gz
            {tabix} -l {out_dir}/sub-set-{num}.sam.tag.gz > {out_dir}/sub-set-{num}.sam.tag.gz.id'''
                     for num, (start, size) in enumerate(break_pos_list)]
        chunk_index = [f'{out_dir}/sub-set-{num}.sam.tag.gz' for num, _ in enumerate(break_pos_list)]

        p.map(os.system, index_cmd, chunksize=1)  # index with multiple threads
        if any(not os.path.exists(f'{i}.csi') for i in chunk_index):
            raise SplitBAMError('Barcode ID blocks not continuous.')

        chunk_tag = []
        for i in chunk_index:
            with open(f'{i}.id') as f:
                chunk_tag.append([i.strip() for i in f])
        chunk_tag = list(dedup(itertools.chain.from_iterable(chunk_tag)))

        final_index = defaultdict(int)
        for key, val in itertools.chain.from_iterable(p.map(get_tag_chunk, chunk_index)):
            final_index[key] += val
        final_index = [(key, final_index[key]) for key in chunk_tag]
        final_index = [(a, b, c) for (a, b), c in zip(final_index, itertools.accumulate(i[-1] for i in final_index))]
        final_index = [(row2[0], row1[-1], row2[1]) for row1, row2 in pairwise(final_index)]  # barcode, start, length
        logging.info(f'Make index: {round(time.monotonic() - t1, 2)} sec.')

        t1 = time.monotonic()
        # split sorted BAM file by tag IDs
        header = f'{out_dir}/sub-set.header'
        os.system(f'{samtools} view -H {out_dir}/sub-set.sam.gz > {header}')
        cmds = [rf'''{bgzip} -b {t1} -s {t2} {out_dir}/sub-set.sam.gz | cat {header} - | \
                    {samtools} view --write-index -o {out_dir}/{cell}.sort.bam''' for cell, t1, t2 in final_index]
        p.map(os.system, cmds, chunksize=10)
        logging.info(f'Split file: {round(time.monotonic() - t1, 2)} sec.')
        os.system(rf'''rm {header} {" ".join(chunk_index)} {" ".join(i + ".csi" for i in chunk_index)}
                    rm {" ".join(i + ".id" for i in chunk_index)}
                    rm {out_dir}/sub-set.sam.gz {out_dir}/sub-set.sam.gz.gzi''')
    logging.info(f'Split {len(chunk_tag) - 1} barcodes ({out_dir}) total: {round(time.monotonic() - start_time, 2)}s.')

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
    split_bam_by_tag(args.bam, args.tag_list, args.out_dir, args.nt, args.tag, args.tag_type)


if __name__ == '__main__':
    set_start_method('spawn')
    main()
