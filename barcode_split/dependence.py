import os
from .exception import SplitBAMError

samtools = '/data/home/lingw/bin/bin/samtools'  # support "view --tag-file" and "sort -t"
bgzip = '/data/home/lingw/bin/bin/bgzip'
mawk = '/data/home/lingw/anaconda3/bin/mawk'
tabix = '/data/home/lingw/bin/bin/tabix'

tools = [samtools, bgzip, mawk, tabix]

for i in tools:
    if not os.path.exists(i):
        script_path = os.path.abspath(__file__)
        raise SplitBAMError(f'{i} not found. Please modify absolute path in {script_path}')

# check "samtools view --tag-file" and "samtools sort -t" 
with os.popen(f'{samtools} view') as h1, os.popen(f'{samtools} sort') as h2:
    m1 = (line.split() for line in h1)
    if not any('--tag-file' in i for i in m1):
        raise SplitBAMError(f'{samtools} view --tag-file not supported.')
    
    m2 = (line.split() for line in h2)
    if not any('-t' in i for i in m2):
        raise SplitBAMError(f'{samtools} sort -t not supported.')
