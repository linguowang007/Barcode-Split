import os
from .exception import SplitBAMError


def get_latest_version_in_path(executable_name):
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        executable_path = os.path.join(path, executable_name)
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            return executable_path
    raise SplitBAMError(f"{executable_name} not found in PATH.")


# get the latest version tools in $PATH
samtools = get_latest_version_in_path("samtools")  # must be supported the "view --tag-file" and "sort -t" options
tabix = get_latest_version_in_path("tabix")
mawk = get_latest_version_in_path("mawk")
bgzip = get_latest_version_in_path("bgzip")
