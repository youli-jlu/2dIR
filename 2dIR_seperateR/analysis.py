import pandas as pd
import numpy as np

freq = pd.read_csv("energy-tot",delimiter=r'\s+')
print(freq.head(5))

print(freq.freq1.mean())
print(freq.freq2.mean())

freq.freq1=freq.freq1-freq.freq1.mean()
freq.freq2=freq.freq2-freq.freq2.mean()

print(freq.head(5))

