try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    from subprocess import check_output
except ImportError:
    import subprocess
    def check_output(*args, **kwargs):
        kwargs['stdout'] = subprocess.PIPE
        p = Popen(*args, **kwargs)
        return p.communicate()[0]
