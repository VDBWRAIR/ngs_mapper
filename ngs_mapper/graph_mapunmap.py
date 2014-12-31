#!/usr/bin/env python

import argparse
import sys
import re
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path

from graph_qualdepth import plot_mapunmap

def parse_args( args=sys.argv ):
    parser = argparse.ArgumentParser(
        description='Graphs Mapped/Unmapped reads for all given qualdepth.json files'
    )

    parser.add_argument(
        dest='jsons',
        nargs='+',
        help='List of qualdepth.json files'
    )

    parser.add_argument(
        '-o',
        dest='outfile',
        default='mapunmap.png',
        help='Filename for the image to be saved'
    )

    return parser.parse_args( args[1:] )

def sample_from_filename( filename ):
    '''
        Returns everything before the .bam in a filename which would hopefully
        be the samplename

        @param filename - Bam filename

        @returns the portion of the filename before .bam
    '''
    m = re.match( '(\S+)\.bam.*', os.path.basename(filename) )
    if not m:
        raise ValueError( "Filename given({}) does not contain .bam".format(filename) )
    return m.group(1)

def get_mapunmap( jsons ):
    '''
        Gets a list of tuples representing mapped,unmapped reads for each sample
        
        @param jsons - List of json files
    
        @returns np.array([(samplename,mapped,unmapped)])
    '''
    samples = []
    mapped_reads = []
    unmapped_reads = []
    for jfile in jsons:
        samples.append( sample_from_filename( jfile ) )
        j = json.load( open(jfile) )
        # Add up all the mapped_reads for every reference
        mreads = 0
        for ref in j:
            if ref != 'unmapped_reads':
                mreads += int(j[ref]['mapped_reads'])
        mapped_reads.append( mreads )
        unmapped_reads.append( int(j['unmapped_reads']) )

    return samples, mapped_reads, unmapped_reads

def plot_samples( ax, smu, mapped_color='g', unmapped_color='r' ):
    '''
        Plots Mapped/Unmapped stacked bar graphs for each sample in smu

        @param ax - Axes to plot on
        @param smu - Sample,MappedReads,UnmappedReads np.ndarray
    '''
    samples = smu[0]
    mapped = smu[1]
    unmapped = smu[2]
    mediani = len(samples)/2.0

    avgm = float(sum(mapped))/len(mapped)
    avgu = float(sum(unmapped))/len(unmapped)

    if isinstance( mediani, float ):
        medianm = mapped[int(mediani)]
        medianu = unmapped[int(mediani)]
    else:
        l = int(mediani)
        u = int(mediani+1)
        medianm = (mapped[l] + mapped[u])/2
        medianu = (unmapped[l] + unmapped[u])/2

    n = np.arange(len(samples))
    mr = ax.bar( n, mapped, color=mapped_color, bottom=unmapped, alpha=0.5 )
    ur = ax.bar( n, unmapped, color=unmapped_color, alpha=0.5 )
    width = mr[0].get_width()
    ax.set_xlim([0,len(samples)])
    ax.set_ylabel( '# Reads' )
    ax.set_xlabel( 'Sample Name' )
    ax.set_xticks( n + width / 2 )
    ax.set_xticklabels( samples, rotation='vertical', fontsize='6', ha='left' )

    # Plot line for the averages of both
    ax.annotate( 'Mapped Reads Avg', xy=(1,avgm) )
    ax.annotate( 'Unmapped Reads Avg', xy=(1,avgu) )
    ax.axhline( y=avgm, color=mapped_color )
    ax.axhline( y=avgu, color=unmapped_color )
    
    # Plot line for the medians of both
    ax.annotate( 'Mapped Reads Median', xy=(mediani,medianm) )
    ax.annotate( 'Unmapped Reads Median', xy=(mediani,medianu) )
    ax.axhline( y=medianm, color=mapped_color, alpha=0.5 )
    ax.axhline( y=medianu, color=mapped_color, alpha=0.5 )

def make_graphic( jsons, outfile ):
    '''
        Creates the graphic and saves it as outfile

        @param jsons - List of json files
        @param outfile - Where to save image
    '''
    smu = get_mapunmap( jsons )

    fig = plt.figure()
    fig.set_size_inches( 20.0, 8.0 )
    gs = gridspec.GridSpec( 1, 2, width_ratios=[20,1] )
    ax1 = plt.subplot( gs[0] )
    ax2 = plt.subplot( gs[1] )

    plot_samples( ax1, smu )
    plot_mapunmap( ax2, sum(smu[1]), sum(smu[2]) )

    fig.savefig( outfile, bbox_inches='tight', dpi=100, pad_inches=0.1 )

def main( args ):
    make_graphic( args.jsons, args.outfile )
