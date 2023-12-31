# Description
The Barcode-Split algorithm is encapsulated within the Python package **barcode_split**.

This algorithm is designed to split BAM files by list of tag values. 

Barcode-split constructs an index of tag values within the SAM file, enabling rapid retrieval of barcode-tagged alignments.
# Dependencies for the barcode_split Package
To use the `barcode_split` package, please ensure the following dependencies are met and that the path of these executable tools are added to the `$PATH` environment variable:
1. python (version 3.6 or higher)
2. samtools (must support the "**samtools view --tag-file**" and "**samtools sort -t**" parameters)
3. bgzip
4. mawk
5. tabix


# Installation of the barcode_split Package

1. Download the source files using "**```git clone https://github.com/linguowang007/Barcode-Split.git```**".
2. Navigate to the **Barcode-Split** folder.
3. Execute "**```python setup.py install```**" in the terminal to complete the installation of the **barcode_split** package.


# How to Use barcode_split

The Barcode-Split algorithm is encapsulated in the `split_bam_by_tag` or `split_bam_by_tag_fast` function within the `barcode_split` package. Both functions have the same usage in terms of parameters and return results. However, the `split_bam_by_tag_fast` function optimizes the sorting procedures using samtools and can leverage multi-core computing capabilities, resulting in faster execution compared to `split_bam_by_tag` in multi-processing models.

Both functions can be utilized in two ways:

## As a Standalone Command-Line Tool. 
After installing **barcode_split** package, the **split_bam_by_tag** tool is automatically added to your system's environment variables. Use "**```which split_bam_by_tag```**" can locate it.

### Usage:
split_bam_by_tag [-h] --bam BAM --tag_list TAG_LIST --out_dir OUT_DIR [--nt NT] [--tag TAG] [--tag_type TAG_TYPE]

split_bam_by_tag_fast [-h] --bam BAM --tag_list TAG_LIST --out_dir OUT_DIR [--nt NT] [--tag TAG] [--tag_type TAG_TYPE]

### Example command in the terminal:

```split_bam_by_tag --bam bamfile.bam --tag_list tags.txt --out_dir output-dir-name --nt 8 --tag CB --tag_type Z```

This command will split bamfile.bam by the values of CB and output the sub-BAM files to output-dir-name.

## As a Python Package Called from Other Python Scripts. 
Add the following code to your Python script (or an interactive Python session):
```
from barcode_split import split_bam_by_tag, split_bam_by_tag_fast

help(split_bam_by_tag)  # help information of split_bam_by_tag function

split_result = split_bam_by_tag(bam="path/to/bamfile.bam",
                                tag_list="path/to/tags.txt",
                                out_dir="out-dir-name",
                                nt=8,
                                tag="CB",
                                tag_type="Z")

# use up to 112 threads, with split_bam_by_tag_fast function 
split_result2 = split_bam_by_tag_fast(bam="path/to/bamfile.bam",
                                      tag_list="path/to/tags.txt",
                                      out_dir="out-dir-name-fast",
                                      nt=112,
                                      tag="CB",
                                      tag_type="Z")
```

# Run test on Demuxlet paper data

```
wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam

samtools index -@ 8 A.merged.bam

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz

gunzip GSM2560245_barcodes.tsv.gz

split_bam_by_tag --bam A.merged.bam --tag_list GSM2560245_barcodes.tsv --out_dir split-test
```
