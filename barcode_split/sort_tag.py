import sys
import argparse
from collections import Counter


def get_cb(*, row, prefix, prefix_len):
    items = row.split()[11:]
    for item in items:
        if item.startswith(prefix):
            return item[prefix_len:], row, len(row)
    return 'no-tag', row, len(row)


def sort_tag(*, out_name, tag='CB', tag_type='Z'):
    prefix = f'{tag}:{tag_type}:'
    prefix_len = len(prefix)
    
    m1 = sorted((get_cb(row=line, prefix=prefix, prefix_len=prefix_len) for line in sys.stdin),
                key=lambda x: (x[0] != 'no-tag', x[0]))
    counter = Counter()
    for cb, line, len1 in m1:
        print(line, end='')
        counter[cb] += len1
    
    with open(out_name, 'wt') as f:
        for barcode in sorted(counter, key=lambda x: (x != 'no-tag', x)):
            print(barcode, counter[barcode], sep='\t', file=f)


# Define a function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Sort BAM lines by tag name. Input from stdin, and output to stdout.')
    parser.add_argument('--out_name', required=True, help='Tag index file name.')
    parser.add_argument('--tag', default='CB', help="Tag name (default: 'CB')")
    parser.add_argument('--tag_type', default='Z', help="Tag values type (default: 'Z')")
    return parser.parse_args()


def main():
    args = parse_arguments()
    # def sort_tag(*, out_name, tag='CB', tag_type='Z'):
    sort_tag(out_name=args.out_name, tag=args.tag, tag_type=args.tag_type)


if __name__ == '__main__':
    main()
