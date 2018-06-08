import pandas as pd
import numpy as np
import glob
import os

class ddict(dict):
    def __missing__(self, key):
        return key

# Read in the MNEc Annotations
# Note: this file can be found here: https://github.com/schae234/PonyTools/blob/master/ponytools/data/MNEc2M_Annotation.csv.gz
info = pd.read_csv("MNEc2M_Annotation.csv.gz",sep=',')
info['legacy'] = np.logical_not([x.endswith('PC') for x in info.MNEcID])

# add an id column so we can do the matching
info['id'] = [f'{r[5]}.{r[6]}_1' for r in info.itertuples()]
info2 = info.copy()
info2['id'] = info2.id.apply(lambda x: x.replace('_1','_2'))
info = pd.concat([info,info2])

# There is a ponytools bug that adds this index column, delete it
del info['Unnamed: 0']  

# Create a dict to hold the results
blasts = dict()
# Get the ensemble mapping to chrids
ensemble_map = ddict()
with open('data/ensemble_id_map.txt') as IN:
    for line in IN:
        k,v = line.strip().split(',')
        ensemble_map[k] = v

# iterate over the blast results in the data dir
for blastresult in glob.glob('MNEc2M_30mer.blast.*'): 
    # Read in the table
    b = pd.read_table(blastresult,names=['id','chrom','perc','start','end'])
    # The SNP position is 35 bases off the end of the blast position
    #b['map_chrom'] = [ensemble_map[x] for x in b.chrom.values]
    #b['map_pos'] = b.end - 35

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
info.to_csv('MNEc2M_30mer.probe_blast_counts.csv')


# to compare number of probe blast counts by genome build
#df = pd.read_csv('MNEc2M_30mer.probe_blast_counts.csv')
#del df['Unnamed: 0']
#n = pd.DataFrame(df.groupby('EquCab2').apply(len),columns=['EC2'])
#n2 = pd.DataFrame(df.groupby('EquCab3').apply(len),columns=['EC3'])
#n['num'] = n.index
#n2['num'] = n2.index
#n = n.merge(n2.reset_index(),left_on='num',right_on='num')
#n['EquCab2'] = n['EC2']
#del n['EC2']
#n['EquCab3'] = n['EC3']
#del n['EC3']

# to get number of variants from which any probe maps uniquely to exactly one location
#EC3_unique = df[df['EquCab3']==1]
#n=pd.DataFrame(EC3_unique.groupby('MNEcID').apply(len),columns=['num_probes'])

# to find the mapping locations for probes
#blast = pd.read_table('MNEc2M_30mer.blast.EquCab3.txt',names=['id','chrom','perc','start','end'])
#blast['pos'] = np.where(blast.id.apply(lambda x: x.split('_')[1]) == '1', blast.start.apply(lambda x: x-1),blast.end.apply(lambda x: x+1))
#del blast['perc']
#EC3_unique = EC3_unique.merge(blast.reset_index(),left_on='id',right_on='id')
#del EC3_unique['index']
#del EC3_unique['ProbeSetID']
#del EC3_unique['AlleleA']
#del EC3_unique['AlleleB']
#del EC3_unique['start']
#del EC3_unique['end']
#del EC3_unique['EquCab2']
#del EC3_unique['EquCab3']
#EC3_unique = EC3_unique.drop_duplicates()

# get SNPs for which both probes map to the same place
#EC3_unique_copy = EC3_unique.copy()
#EC3_unique_copy = EC3_unique_copy[EC3_unique_copy.duplicated()]
#EC3_unique_copy['EquCab2_chrom'] = EC3_unique_copy['chrom_x']
#del EC3_unique_copy['chrom_x']
#EC3_unique_copy['EquCab2_pos'] = EC3_unique_copy['pos_x']
#del EC3_unique_copy['pos_x']
#EC3_unique_copy['EquCab3_chrom'] = EC3_unique_copy['chrom_y']
#del EC3_unique_copy['chrom_y']
#EC3_unique_copy['EquCab3_pos'] = EC3_unique_copy['pos_y']
#del EC3_unique_copy['pos_y']
#EC3_unique_copy.to_csv('MNEc2M_30mer.both_probes_agree.csv')
# just the ones mapped to the same chromosome as before
#EC3_unique_copy[EC3_unique_copy['EquCab2_chrom']==EC3_unique_copy['EquCab3_chrom']].to_csv('MNEc2M_30mer.both_probes_agree_same_chrom.csv')

# get the SNPs mapping uniquely to exactly one place
#del EC3_unique['id']
#EC3_unique = EC3_unique.drop_duplicates()
#n = pd.DataFrame(EC3_unique.groupby('MNEcID').apply(len),columns=['num'])
#n = n[n['num']==1]
#n = n.reset_index()
#EC3_unique=EC3_unique.merge(n.reset_index(),left_on='MNEcID',right_on='MNEcID')
#del EC3_unique['index']
#del EC3_unique['num']
#EC3_unique.to_csv('MNEc2M_30mer.unique_remap.csv')
