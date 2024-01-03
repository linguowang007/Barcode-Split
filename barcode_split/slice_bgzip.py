import os
import bisect
import argparse
from array import array
from .dependence import bgzip


def slice_bgzip(*, bgzip_file, pos1, pos2, out_file, file_mode='ab'):
    """
    Extracts a slice from a bgzip compressed file.

    Args:
       bgzip_file (str): Path to the bgzip compressed file.
       pos1 (int): Starting byte position of the slice.
       pos2 (int): Ending byte position of the slice.
       out_file (str): Output file path.
       file_mode (str): File mode for output file. Either 'ab' (append binary) or 'wb' (write binary).

    Returns:
    str: The path to the output file.

    Raises:
       ValueError: If file_mode is not 'ab' or 'wb'.
    """
    if file_mode not in ['ab', 'wb']:
        raise ValueError("file_mode must be 'ab' or 'wb'")

    index = f'{bgzip_file}.gzi'
    with open(index, 'rb') as ori_f:
        data = array('q')
        data.fromfile(ori_f, 1)
        data.fromfile(ori_f, data[0] * 2)

    gzip_index = data[1::2]
    txt_index = data[2::2]

    insert1 = bisect.bisect_left(txt_index, pos1)
    insert2 = bisect.bisect_left(txt_index, pos2)

    con_gzip = (gzip_index[insert1], gzip_index[insert2 - 1],)
    left_pos = (pos1, txt_index[insert1])
    right_pos = (txt_index[insert2 - 1], pos2)

    a, b = left_pos
    if file_mode == 'ab':
        os.system(f'{bgzip} -b {a} -s {b - a} {bgzip_file} | {bgzip} >> {out_file}')
    else:
        os.system(f'{bgzip} -b {a} -s {b - a} {bgzip_file} | {bgzip} > {out_file}')

    with open(bgzip_file, 'rb') as ori_f:
        block_size = 1024 * 4
        a, b = con_gzip
        ori_f.seek(a)
        q_len = b - a
        num, len1 = divmod(q_len, block_size)

        count = 0
        with open(out_file, 'ab') as write_f:
            while count < num:
                con = ori_f.read(block_size)
                count += 1
                write_f.write(con)
            con = ori_f.read(len1)
            write_f.write(con)

    a, b = right_pos
    os.system(f'{bgzip} -b {a} -s {b - a} {bgzip_file} | {bgzip} >> {out_file}')

    return out_file


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract a slice from a bgzip compressed file.")
    parser.add_argument("--bgzip_file", required=True, help="Path to the bgzip compressed file.")
    parser.add_argument("--pos1", required=True, type=int, help="Starting byte position of the slice.")
    parser.add_argument("--pos2", required=True, type=int, help="Ending byte position of the slice.")
    parser.add_argument("--out_file", required=True, help="Output file path.")
    parser.add_argument("--file_mode", default='ab', choices=['ab', 'wb'],
                        help="File mode for writing to output file (default: 'ab').")
    return parser.parse_args()


def main():
    args = parse_arguments()
    slice_bgzip(bgzip_file=args.bgzip_file, pos1=args.pos1, pos2=args.pos2, out_file=args.out_file,
                file_mode=args.file_mode)


if __name__ == "__main__":
    main()
