#!/usr/bin/env python

import json
import argparse
import sys
from os.path import basename
import numpy as np

def main( args ):
    if args.title is None:
        title = basename(args.outfile)
    else:
        title = args.title
    make_graphic( args.jsonfile, args.outfile, args.ref, titleprefix=title )

def plot_depths( ax, xvals, yvals, maxdepth, color, title ):
    ax.set_title( "{0} Depth/Qual".format(title) )
    ax.set_xlabel( "Reference Position" )
    ax.set_ylabel( "Depth" )
    ax.fill_between( xvals, yvals, facecolor=color, alpha=0.5 )
    #ax.plot( xvals, yvals, c=color )
    ax.set_xlim([0,xvals[-1]])
    ax.set_ylim([0,maxdepth])
    ax.minorticks_on()

def plot_quals( ax, xvals, yvals, color ):
    from matplotlib.ticker import MultipleLocator
    from matplotlib.ticker import LinearLocator
    ticks = None
    # It looks nicer if you have ticks on multiples of 500,
    # but if the genome is very large just do every 20
    if len( xvals ) < 15000:
        # Put ticks every 500 on the x axis
        ticks = MultipleLocator( 500 )
    else:
        ticks = LinearLocator( 10 )
    ax.xaxis.set_major_locator( ticks )
    #ax.plot( xvals, yvals, c=color )
    ax.fill_between( xvals, yvals, facecolor=color, alpha=0.2 )
    ax.set_ylabel( "Quality" )
    ax.set_xlim([0,xvals[-1]])
    ax.set_ylim([0,40])
    ax.minorticks_on()
    ax.axhline( y=37, color='g', linestyle=':' )

def plot_mapunmap( ax, mapped, unmapped, mapped_color='g', unmapped_color='r', mapped_read_text_color='black', unmapped_read_text_color='white', ylabel='Unmapped/Mapped Reads'):
    mr = ax.bar( [0], [mapped], 0.05, color=mapped_color, bottom=unmapped )
    ur = ax.bar( [0], [unmapped], 0.05, color=unmapped_color )
    ax.minorticks_off()
    ax.set_autoscalex_on(False)
    ax.set_xlim([0,0.05])
    ax.set_ylim([0,mapped+unmapped])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel( ylabel )

    def autolabel( rects, color=mapped_read_text_color, bottom=0, text_under=False ):
        # Text size in pt
        text_size = 12.0
        # DPI of the figure
        dpi = ax.figure.dpi
        # The text size in pixels
        text_size_px = text_size * 72 / dpi
        for rect in rects:
            height = rect.get_height()
            text_left = rect.get_x()+rect.get_width()/2.0
            if text_under:
                scale = 1 - text_size*2.5 / 1000
            else:
                scale = 1 + text_size / 1000
            text_bottom = bottom + height * scale
            ax.text( text_left, text_bottom, '%d'%int(height), size=text_size, ha='center', va='bottom', color=color )

    autolabel( mr, bottom=unmapped )
    autolabel( ur, text_under=True, color=unmapped_read_text_color )

def make_graphic( qualdepthfile, outputfile, ref=None, titleprefix='', compress_data=50 ):
    '''
        Makes a graphic for a reference showing depth and avg qualities

        @param qualdepthfile - Should be qualdepth.json file
        @param outputfile - Where to save the image
        @param ref - Which reference to do the image for
        @param titleprefix - What to put in title before the Qual/Depth text
        @param compress_data - How many ref positions to skip between each data point
        
        NOTE: This works for single reference mapped only!!
    '''
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np

    # Load the json
    j = json.load( open(qualdepthfile) )

    refs = [r for r in j.keys() if r != 'unmapped_reads']
    if ref is None:
        ref = refs[0]

    # Find the maximum depth so that we can standardize the yaxis
    max_depth = 0
    for r in j:
        if r != 'unmapped_reads':
            maxd = j[r]['maxd']
            if maxd > max_depth:
                max_depth = maxd

    fig = plt.figure()
    fig.set_size_inches( 20.0, 8.0 )
    gs = gridspec.GridSpec(1,2, width_ratios=[20,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = ax1.twinx()

    depths = j[ref]['depths'][::compress_data]
    quals = j[ref]['avgquals'][::compress_data]
    mapped_reads = j[ref]['mapped_reads']
    unmapped_reads = j['unmapped_reads']
    xvals = range(0,len(j[ref]['depths']), compress_data)

    plot_depths( ax1, xvals, depths, max_depth, title=titleprefix, color='blue' )
    plot_quals( ax3, xvals, quals, color='green' )
    plot_mapunmap( ax2, mapped_reads, unmapped_reads )
    #TODO: Here, get and plot paired v. unpaired
#    paired, unpaired = samtools.get_pair_info()
#    plot_mapunmap( ax2, paired, unpaired, mapped_color='blue', unmapped_color='r', mapped_read_text_color='black', unmapped_read_text_color='white', ylabel='Paired/Unpaired Reads'):

    fig.savefig( outputfile, bbox_inches='tight', dpi=100, pad_inches=0.1 )

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Make a Depth/Quality graphic from a qualdepth.json file from running bam_to_qualdepth on a bam file'
    )

    parser.add_argument(
        'jsonfile',
        help='Output file from running bam_to_qualdepth'
    )

    parser.add_argument(
        '-r',
        '--reference',
        dest='ref',
        default=None,
        help='What reference to generate the image for[Default: First found reference]'
    )

    default_output='qualdepth.png'
    parser.add_argument(
        '-o',
        '--outfile',
        dest='outfile',
        default=default_output,
        help='Where to save the image file[Default: %(default)s]'
    )

    parser.add_argument(
        '-t',
        '--title',
        dest='title',
        default=None,
        help='Title of the graphic[Default: basename of outfile argument]'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
