import shutil
from os.path import *
import os
from glob import glob

import log
logger = log.setup_logger( __name__, log.get_config() )

def sync_sanger( runpath, ngsdata ):
    pass

def sync_run( runpath, ngsdata ):
    rawd = join( ngsdata, 'RawData', 'Sanger' )
    run = basename(runpath)
    rund = join( rawd, run )

    if not isdir( rund ):
        os.makedirs( rund )

    reads_to_copy = glob( join( runpath, '*.ab1' ) )
    print reads_to_copy
    for read in reads_to_copy:
        rp = join( rund, basename(read) )
        if not exists( rp ):
            logger.info( 'Copied {} to {}'.format(read,rp) )
            shutil.copy( read, rp )
        else:
            logger.info( '{} already existed so it was skipped'.format(rp) )
    
def sync_readdata( rawdir, ngsdata ):
    pass



















'''
 for r in reads:
...   sn = re.match( '(\S+?)_[FR]\d*', basename(r) ).groups(1)[0]
...   run = dirname( r )
...   lnk = relpath( r, join('/home','EIDRUdata','NGSData', 'ReadsBySample', sn ) )
...   sd = join( '/home', 'EIDRUdata', 'NGSData', 'ReadsBySample', sn )
...   if not isdir( sd ):
...     os.makedirs( sd )
...   os.symlink( lnk, join( sd, basename(r) ) )

'''
