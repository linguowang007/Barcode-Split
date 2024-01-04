from setuptools import setup, find_packages

setup(
    name='barcode_split',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'split_bam_by_tag = barcode_split.split_bam_file:main',
            'slice_bgzip = barcode_split.slice_bgzip:main',
            'split_bam_by_tag_fast = barcode_split.split_bam_by_tag_fast:main'
        ]
    },
)
