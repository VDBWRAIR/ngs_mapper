import os
import os.path

def build_datafiles(prefix, datadir):
    '''
    Return a data_files like structure for distutils for a given 
    data dir that

    datadir is relative to setup.py

    :param str prefix: Where to install into
    :param str datadir: The datadir to iterate through
    :return: [(prefix, [files]), (prefix/dir1, [dir1's files])]
    '''
    manifest = []
    for root, dirs, files in os.walk(datadir):
        _prefix = os.path.normpath(root.replace(datadir, prefix))
        manifest.append((_prefix, [os.path.join(root,f) for f in files]))
    return manifest
