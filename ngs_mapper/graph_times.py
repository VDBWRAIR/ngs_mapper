from glob import glob
from os.path import *
import numpy as np
import matplotlib.pyplot as plt

from nose.tools import ok_, eq_

from datetime import datetime
import log

logc = log.get_config()
logger = log.setup_logger( 'graph_times', logc )

def main():
    ss = start_stop( 'Projects' )
    logger.info( "Plotting all projects inside of {0}".format('Projects') )
    x,y = [],[]
    samplenames = sorted(ss.keys())
    for sn in samplenames:
        x.append( sn )
        y.append( ss[sn] )
    fig = plt.figure()                                                                                                                                                                                                                                                                
    fig.set_size_inches( 20.0, 8.0 )
    fig.suptitle( 'Pipeline Time per Sample' )
    ax = plt.gca()
    ax.plot( range(len(x)), y )
    ax.set_xlim([0,len(x)-1])
    ax.set_ylim([0,max(y)])
    ax.set_xticks( range(0,len(x)) )
    ax.set_xticklabels( x, rotation='vertical' )
    ax.set_ylabel( 'Seconds' )
    plt.savefig( 'PipelineTimes.png', bbox_inches='tight', dpi=100, pad_inches=0.1 )

def datediff( start_stop ):
    '''
        >>> from datetime import timedelta
        >>> start = '2014-03-18 14:51:41,000'
        >>> stop = '2014-03-18 14:59:26,000'
        >>> ss = (start,stop)
        >>> eq_( 465 , datediff( ss ) )
    '''
    fmt = '%Y-%m-%d %H:%M:%S,%f'
    start = datetime.strptime( start_stop[0], fmt )
    stop = datetime.strptime( start_stop[1], fmt )
    return (stop - start).seconds

def start_stop( basedir ):
    '''
        >>> ss = start_stop( 'Projects' )
        >>> k = ss.keys()[0]
        >>> ok_( isinstance( ss[k], int ) )
    '''
    projects = get_projects( basedir )
    ss = {}
    for p in projects:
        proj = basename( p )
        diff = datediff( start_stop_for_project( p ) )
        if diff < 60:
            logger.warning( "{0} ran in only {1} seconds".format(p,diff) )
        ss[proj] = diff
    return ss

def start_stop_for_project( projectpath ):
    '''
        >>> p = get_projects( 'Projects' )[0]
        >>> ss = start_stop_for_project( p )
        >>> eq_( (1,1), ss )
    '''
    import re
    p = '(\d{4}-\d{2}-\d{2}\s\d{2}:\d{2}:\d{2},\d+).*?--- (Starting|Finished)'
    log = join( projectpath, basename( projectpath ) + '.log' )
    ss = None
    with open( log ) as fh:
        ss = re.findall( p, fh.read() )
    return ss[0][0], ss[1][0]

def get_projects( basedir ):
    '''
        >>> pds = get_projects( 'Projects' )
        >>> eq_( len( pds ), 94 )
    '''
    return glob( join( basedir, '*' ) )
