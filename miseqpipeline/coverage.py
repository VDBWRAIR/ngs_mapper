"""
Generates a single image that has a sub-image for every reference found in all of the projects supplied.
It uses the already generated :doc:`qualdepthjson` for each sample to create a single multi-colored line to easily show where Gaps, LowCoverage, LowQuality and LowCoverage+LowQuality regions are.

Help usage
==========

    .. code-block:: bash

        sample_coverage --help

Examples
--------

Typically the pipeline is run such that there is a Projects directory that contains the analysis and all of the files for each sample that was run.
You can easily manually run sample_coverage on any number of project directories as follows:

Run specific projects
---------------------

    .. code-block:: bash

        sample_coverage Projects/sample1 Projects/sample2 ...

Run on all project directories
------------------------------

    .. code-block:: bash

        sample_coverage Projects/*

Run for specific references
---------------------------

Sometimes you may not want all references compiled into the resulting image(more than likely since the image will become very blury as the script will have to make it not as good quality)

Exclude list
------------

You can specify an exclude list to exclude keywords that are in your references as follows

Exclude, pH1N1 and H3N2 references:

    .. code-block:: bash

        sample_coverage Projects/* --exclude pH1N1 H3N2

Include list
------------

You can alternatively specify a list of names only to include which has precedence over the exclude list
Exclude pH1N1 H3N2, but show all references with /MP/ in them(even pH1N1 and H3N2)

    .. code-block:: bash

        sample_coverage Projects/* --exclude pH1N1 H3N2 --include '/MP/'
"""
from glob import glob
import json
from os.path import join, basename
from collections import defaultdict, OrderedDict
import argparse

from bqd import (
    lines2d_from_regions,
    regions_from_qualdepth,
    G, N, LQ, LC, LCQ,
    REGIONTYPES,
)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

# Line color mappings for region types
LINEARGS = {
    G: {'color':'red', 'linewidth':5},
    N: {'color':'grey', 'linewidth':5},
    LQ: {'color':'blue', 'linewidth':5},
    LC: {'color':'orange', 'linewidth':5},
    LCQ: {'color':'yellow', 'linewidth':5},
}

def filter_refs(refs, includes, excludes):
    '''
    Filter references using includes and excludes string lists
    Includes take precedence over excludes
    '''
    # Here will be are unique references
    allrefs = set()
    # Loop through references
    for ref in refs:
        # only include what was requested
        if includes:
            for include in includes:
                if include in ref:
                    allrefs.add(ref)
        else:
            for exclude in excludes:
                if exclude not in ref:
                    allrefs.add(ref)
        if not includes and not excludes:
            # No filters so do all
            allrefs.add(ref)
    return allrefs

def refs_from_project(projectpath, includes, excludes):
    '''
    Get a set of uniq references from a project path

    includes and excludes get passed to filter_refs to filter
     the refs down if wanted

    returns a set object
    '''
    # The uniq references after filtering
    refs = set()
    qualdepths = load_project_qualdepth(projectpath)
    # Add filtered references for this qualdepth
    return filter_refs(qualdepths.keys(), includes, excludes)

def get_allrefs(projects, includes, excludes):
    '''
    Get a unique set of reference names from all qualdepth files in all project
    directories
    '''
    # Start with empty set then fill it in
    allrefs = set()
    # Loop all qualdepth.json and get all reference names
    # Some projects may be mapped to diff refs
    for p in projects:
        # Update master set
        allrefs.update(refs_from_project(p, includes, excludes))
    return allrefs

def load_project_qualdepth(projpath):
    '''
    Simply load the qualdepth.json file for a given project path
    '''
    try:
        qualdepthfile = glob(join(projpath, '*.bam.qualdepth.json'))[0]
    except IndexError as e:
        raise ValueError('{0} missing qualdepth file'.format(projpath))
    return json.load(open(qualdepthfile))

def get_perreference_from_projects(projects, allrefs, refax, gap, lowqual, lowcov, lineargs):
    '''
    Get a dictionary keyed by each reference in allrefs
    Each item contains a list of generators that generate Line2D objects
    for a sample's regions

    projects - list of project paths
    allrefs - sequence of reference names
    refax - {ref1:axes, ref2:axes,...}
    gap, lowqual, lowcov are integers dictating how to call coverage regions
    lineargs is a dictionary of kwargs for Line2D
    '''
    # Now build Line2D for everything
    # This builds lines for each sample in each reference
    perreference = defaultdict(list)
    for sampleno, projdir in enumerate(projects, start=1):
        qualdepths = load_project_qualdepth(projdir)
        # Add each sample's region lines to the reference key
        for ref in allrefs:
            if ref not in qualdepths:
                # Skip missing references for this sample
                continue
            # Get correct plot
            ax = refax[ref]
            # Get the qualdepth for reference
            qualdepth = qualdepths[ref]
            ax.set_xlim(0,int(qualdepth['length']))
            regions = regions_from_qualdepth(qualdepth, gap, lowqual, lowcov)
            # Each reference will hold a list of line2d generators
            perreference[ref].append(lines2d_from_regions(sampleno, regions, lineargs))
    return perreference

def set_figure_size(perreference, figure):
    '''
    Set figure size and dpi based on number of references and number of
    samples
    '''
    # Get maximum number of samples in any reference
    maxnumsamples = max([len(perreference[ref]) for ref in perreference])
    # Set the size of the figure in inches
    # Seems that 0.375 inches per sample is enough
    # Max num samples per referece * num references/2 * 0.357
    # num references / 2 because 2 columns
    numrefs = len(perreference)
    length = maxnumsamples * (numrefs / 2.0) * 0.375
    if length > 400:
        # Scale dpi down
        scale = 400.0 / length
        dpi = figure.get_dpi()
        dpi = scale * dpi
        print "Large amount of samples, scaling image quality down"
        figure.set_dpi(dpi)
        # Reset length to our 400 maximum
        length = 400
        bbox = figure.get_window_extent().transformed(figure.dpi_scale_trans.inverted())

    if length < 1:
        length = 2
    figure.set_size_inches(20.0, length)

def plot_all_references(projects, perreference, refax):
    # Samplenames
    # Should be just all the basenames of projects
    # 0 index should be blank because the enumerate is start=1 above
    samplenames = [''] + [basename(p) for p in projects]

    # Do all plotting now
    for ref in perreference:
        # Get axes
        ax = refax[ref]
        # Set the title of the axes to current reference
        ax.set_title(ref)
        # Set the x-axis label
        ax.set_xlabel('Reference Position')

        # Plot the reference
        samples = plot_reference(ref, perreference[ref], samplenames, ax)

        # Set Y-Axis details
        # Set y min,max
        ax.set_ylim(0,len(samples))
        # Tell y-axis to have same number of ticks as # samples
        ax.set_yticks(xrange(len(samples)))
        # Set each tick label
        ax.set_yticklabels(samples)
        # Set the Label that shows up for Y-Axis
        ax.set_ylabel('Samples')

def plot_reference(ref, line2dseq, samplenames, ax):
    '''
    Plot reference on a given axes using line2dsequences for reference
    '''
    # Samples for this axis in correct order
    samples = ['']
    # Also keep a quicker lookup
    samplesseen = set()
    # perreference contains key: list(generatorofline2ds, generatorofline2ds,)
    for y, sampleslines in enumerate(line2dseq, start=1):
        # Get samplename from samplenames
        # Iterate all the line2d objects for each sample
        for line in sampleslines:
            # First get index which is same as Line2D ydata
            i = line.get_ydata()[0]
            # Fetch sample from samplenames list
            sample = samplenames[i]
            # Add sample to ylabels list(samples)
            if sample not in samplesseen:
                samples.append(sample)
                samplesseen.add(sample)
            # Reset this line to the correct y position
            line.set_ydata([y,y])
            # Add the line to the axes
            ax.add_line(line)
    return samples

def create_legend(figure, regiontypes, lineargs):
    # Draw legend
    # First get a mock regions so we can get line2d objects for each region type
    legendlines = [
        Line2D([0,0],[0,0],**lineargs[regiontype]) for regiontype in regiontypes
    ]
    # This is all magic and wizardry here
    # I don't get how it works, I just try things until
    # it shows the way I want
    figure.legend(
        legendlines,
        regiontypes,
        bbox_to_anchor=(1,1),
        borderaxespad=1,
        loc='upper left',
        bbox_transform=figure.transFigure,
        #ncol=len(legendlines),
    )

def create_figure_for_projects(projects, includes, excludes, lineargs, regionmins):
    gap, lowqual, lowcov = regionmins

    # Will hold our unique reference names from all projects
    allrefs = get_allrefs(projects, includes, excludes)

    # Get the figure object
    fig = plt.figure()

    # Grid of figures for each reference
    # 2 columns, numrefs rows
    gs = gridspec.GridSpec(len(allrefs), 2, width_ratios=[1,1])

    # Map reference name to it's subplot
    # Keep sorted order of references so they are grouped together
    refax = OrderedDict([(ref,plt.subplot(gs[i])) for i,ref in enumerate(sorted(allrefs))])

    # Get line segments for each sample broken down by reference
    perreference = get_perreference_from_projects(projects, allrefs, refax, gap, lowqual, lowcov, lineargs)

    #Setup figure size
    set_figure_size(perreference, fig)

    plot_all_references(projects, perreference, refax)

    create_legend(fig, REGIONTYPES, lineargs)

    # Adjust the subplots(reference graphs) so they have
    # correct spacing
    gs.tight_layout(fig)
    return fig

def parse_args():
    parser = argparse.ArgumentParser(
        'Create graphic to show coverage per sample broken down by reference name'
    )

    parser.add_argument(
        'projects',
        nargs='+',
        help='List of projects to generate from. Each project needs a sample.bam.qualdepth.json'
    )

    parser.add_argument(
        '--output',
        default='Coverage.png',
        help='Path to save output image to[Default: %(default)s]'
    )

    parser.add_argument(
        '--lowcov',
        default=10,
        type=int,
        help='< this is low coverage[Default: %(default)s]'
    )

    parser.add_argument(
        '--lowqual',
        default=25,
        type=int,
        help='< this is low quality[Default: %(default)s]'
    )

    parser.add_argument(
        '--gap',
        default=0,
        type=int,
        help='< this is low quality[Default: %(default)s]'
    )

    parser.add_argument(
        '--exclude',
        nargs='+',
        default=[],
        help='List of strings to match in reference names to exclude ' \
            'from graphic[Default: %(default)s]'
    )

    parser.add_argument(
        '--include',
        nargs='+',
        default=[],
        help='List of strings to match in reference names to include ' \
            'in grpahic[Default: %(default)s]'
    )

    return parser.parse_args()

def main():
    args = parse_args()

    lowcov = args.lowcov
    lowqual = args.lowqual
    gap = args.gap
    lineargs = LINEARGS
    projects = sorted(args.projects)
    excludes = ['unmapped_reads'] + args.exclude
    includes = args.include

    regionmins = [gap,lowqual,lowcov]

    fig = create_figure_for_projects(projects, includes, excludes, lineargs, regionmins)
    try:
        fig.savefig(args.output, dpi=fig.dpi, bbox_inches='tight')
    except ValueError as e:
        print "!!!!!!!!!! Error: Image size too large to create !!!!!!!!!!!!!!"
        raise e

if __name__ == '__main__':
    main()
