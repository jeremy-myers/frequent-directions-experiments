import sys
from utils.syntheticDataMaker import SyntheticDataMaker

n = 1000
d = 100
k = 5

fname = "matrix.txt"
file  = open( fname, "w" )

# this is only needed for generating input vectors
dataMaker = SyntheticDataMaker()
dataMaker.initBeforeMake(d, k, signal_to_noise_ratio=10.0)

for i in xrange(n):
    row = dataMaker.makeRow()
    dataMaker.writeToFile( row, file )

file.close()
