import os
from barcode_split import split_bam_by_tag
from multiprocessing import set_start_method

if __name__ == '__main__':
    set_start_method('spawn')

    cmd = r'''
    wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam
    samtools index -@ 8 A.merged.bam
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz
    gunzip GSM2560245_barcodes.tsv.gz
    '''
    os.system(cmd)

    bam_file = 'A.merged.bam'
    tag_list_file = 'GSM2560245_barcodes.tsv'
    out_dir_name = 'split_test'

    result = split_bam_by_tag(bam_file, tag_list_file, out_dir_name, nt=16)

    print(f'Split {len(result)} cells.')

    for cell, sub_bam in result.items():
        print(cell, sub_bam, sep='\t')
