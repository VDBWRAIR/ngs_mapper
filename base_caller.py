


from stats_at_refpos import stats
import argparse

# arg to enter bam file
parser = argparse.ArgumentParser(description="Bam_file")

parser.add_argument(
	dest='bam',
	help='Bam file path'
)



# Alignment file
alignment_bam = args.bam


# stats( bamfile, regionstr, minmq, minbq, maxd ):

stats( {}, regionstr, 1, 1, 100000).format(alignment_bam):

def label_N(stats, minbq = 25):
	stats2 = {}
	stats2['depth'] = stats['depth']
	stats2['mqualsum'] = stats['mqualsum']
	stats2['bqualsum'] = stats['bqualsum']


	for base, quals in stats.iteritems():
        	# Only interested in base stats in this loop
        	if base not in ('depth','mqualsum','bqualsum'):

			# generates a list called bquals
			bquals = quals['baseq']
			

			# loop to examine the quality score and identifyes bases with a quality score less than the minbq of 25
			for q in bquals:
				k = base
				if q < minbq:
					k = 'N'

				# adds the N to the nucleotides (A C G T and N)
				if k not in stats2:
					stats2[k] = {'baseq':[]}
				stats2[k] ['baseq'].append( q )


	return stats2

# call the nucleotides
def caller( {}, regionstr, 1, 1, 100000, mindp=10, min_th=0.8).format(alignment_bam):

	stats = stats(bamfile, regionstr, minmq, minbq, maxd):
	stast2 = label_N(stats, minbq)


	# if the quality is 25 or greater ignore the N bases
	if stats2['depth'] >= mindp:
		del stats2['N']

	else:
		# defines if the base is an N
		np = len(stats2['N']['baseq'])/(stats2['depth']*1.0)
		if np > (1-min_th)
			return 'N'

	return call_on_pct(stats2, min_th)






			

			
