#!/home/users/jtm1171/.conda/envs/qiime2-amplicon-2024.5/bin/python3.9
import sys
import pandas as pd
from qiime2 import Metadata

# first we import the metdata into pandas
mf = pd.read_csv(sys.argv[1], sep='\t',index_col=0)
# next we aggregate by subjects (i.e. 'host_subject_id') 
# and keep the first instance of 'diagnosis_full' by subject.
mf = mf.groupby('Site_Corrected').agg({'NERR':'first','Site_Corrected':'first','salinity':'first'})
# now we save the metadata in QIIME2 format.
mf.index.name = '#SampleID'
outname = sys.argv[1].split('.')[0] + '-subject-metadata.tsv'
mf.to_csv(outname, sep='\t')
quit()
