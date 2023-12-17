wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam
samtools index -@ 8 A.merged.bam

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz
gunzip GSM2560245_barcodes.tsv.gz

split_bam_by_tag --bam A.merged.bam --tag_list GSM2560245_barcodes.tsv --out_dir split-test --nt 16
