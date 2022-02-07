import sys
from numpy.linalg import norm
from numpy import dot
from numpy import loadtxt
from numpy import savetxt

from utils.syntheticDataMaker import SyntheticDataMaker
from frequentDirections import FrequentDirections

n = 1000
d = 100
ell = 20
k = 5
print_level = 0

matrix = loadtxt( "matrix.txt" )

shrink = True
shrink = False

if shrink == True:
    save_fname = "shrinkt.txt"
else:
    save_fname = "shrinkf.txt"

# This is where the sketching actually happens
sketcher = FrequentDirections(d,ell,shrink,print_level)
for i in xrange(n):
    row = matrix[i,:]
    sketcher.append(row)
sketch = sketcher.get()

savetxt( save_fname, sketch, delimiter=" " )
