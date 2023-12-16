from setuptools import setup, find_packages

setup(
    name='bebam',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'split_bam_by_tag = bebam.split_bam_file:main'
        ]
    },
)
