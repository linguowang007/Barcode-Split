import os
from .exception import SplitBAMError

samtools = rf'/data/home/lingw/bin/bin/samtools'  # version 1.16.1 or higher, support "view --tag-file"
bgzip = rf'/data/home/lingw/bin/bin/bgzip'
mawk = '/data/home/lingw/anaconda3/bin/mawk'
tabix = '/data/home/lingw/bin/bin/tabix'

tools = [samtools, bgzip, mawk, tabix]
for i in tools:
    if not os.path.exists(i):
        raise SplitBAMError(f'{i} not found.')
