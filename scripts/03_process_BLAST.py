import pandas as pd
import numpy as np
import glob
import os

class ddict(dict):
    def __missing__(self, key):
        return key

# Read in the MNEc Annotations
# Note: this file can be found here: https://github.com/schae234/PonyTools/blob/master/ponytools/data/MNEc2M_Annotation.csv.gz
info = pd.read_csv("/home/rob/Codes/PonyTools/ponytools/data/MNEc2M_Annotation.csv.gz",sep=',')
info['legacy'] = np.logical_not([x.endswith('PC') for x in info.MNEcID])

# add an id column so we can do the matching
info['id'] = [f'{r[5]}.{r[6]}' for r in info.itertuples()]

# There is a ponytools bug that adds this index column, delete it
del info['Unnamed: 0']  

# Create a dict to hold the results
blasts = dict()
# Get the ensemble mapping to chrids
ensemble_map = ddict()
with open('../data/ensemble_id_map.txt') as IN:
    for line in IN:
        k,v = line.strip().split(',')
        ensemble_map[k] = v

# iterate over the blast results in the data dir
for blastresult in glob.glob('../data/MNEc2M.blast.*'): 
    # Read in the table
    b = pd.read_table(blastresult,names=['id','chrom','perc','start','end'])
    # The SNP position is 35 bases off the end of the blast position
    b['map_chrom'] = [ensemble_map[x] for x in b.chrom.values]
    b['map_pos'] = b.end - 35

    b['probe_chrom'] = b.id.apply(lambda x: x.split('.')[0])
    b['probe_pos'] = b.id.apply(lambda x: x.split('.')[1])

    del b['start']
    del b['end']
    del b['chrom']

    # extract the ref genome name fro the file
    file_source = os.path.basename(blastresult).split('.')[2]
    # save the blast results in a dictionary for later use
    blasts[file_source] = b

    # Group by ID which will and iterate over results
    # apply the length function to each grouping to get the number of places
    # that probe BLASTED to
    n = pd.DataFrame(b.groupby('id').apply(len),columns=[file_source])
    # Add the results to the info data frame using the merge function
    info = info.merge(n.reset_index(),left_on='id',right_on='id')




# Save the table to a file
#info.to_csv('../data/MNEc2M.probe_blast_counts.csv')
