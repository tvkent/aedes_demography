import sys
import pandas as pd
import numpy as np

"""
0: chromosomes of sweepfinder results
1: outpath
"""

test=pd.read_csv(sys.argv[1],sep='\t')

test2=pd.read_csv(sys.argv[2],sep='\t')

test3=pd.read_csv(sys.argv[3],sep='\t')

test_merged=pd.concat([test,test2,test3],ignore_index=True)

selset=test_merged.sort_values('LR',ascending=False).head(int(len(test_merged['LR'])*0.01))[['chrom','location']]

selset['start']=selset['location']-10001
selset['end']=selset['location']+10000
selbed=selset[['chrom','start','end']]
selbed=selbed.sort_values(['chrom','start'])

selbed.to_csv(sys.argv[4],sep="\t", index=False, header=False)
