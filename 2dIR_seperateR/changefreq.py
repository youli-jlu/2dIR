import pandas as pd
import numpy as np
import itertools as it

rawdata = pd.read_csv("IRgap1.dat",delimiter=r'\s+')
print(rawdata.head(5))

i,j=0,0
for index,raw in rawdata.iterrows():
    if (raw['w2']==2459.26):
        i+=1
        j=0
    raw['w1']=i
    raw['w2']=j
    j+=1

print(rawdata.head(5))
rawdata.to_csv("IRgap1-grid.dat",sep=" ",index=False)


