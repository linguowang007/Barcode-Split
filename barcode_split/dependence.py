import os
from .exception import SplitBAMError

samtools = '/data/home/lingw/bin/bin/samtools'  # support "view --tag-file" and "sort -t"
bgzip = '/data/home/lingw/bin/bin/bgzip'
mawk = '/data/home/lingw/anaconda3/bin/mawk'
tabix = '/data/home/lingw/bin/bin/tabix'

tools = [samtools, bgzip, mawk, tabix]
script_path = os.path.abspath(__file__)

for i in tools:
    if not os.path.exists(i):
        raise SplitBAMError(f'{i} not found. Please modify variable in {script_path}')
