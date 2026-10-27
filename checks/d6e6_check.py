# Table 3: D6 = T(1,1,3) and E6 = T(1,2,2) at the four subcubic cells and the exponential jet (exact direct edge-algebra computation).
# Usage: python3 d6e6_check.py ../certificates
import sys
sys.path.insert(0, sys.argv[1])
from direct_hilbert import hilbert
P=0  # exact rational arithmetic (a prime here would compute modulo it)
D6=[(1,2),(2,3),(3,4),(4,5),(2,6)]
E6=[(1,2),(2,3),(3,4),(4,5),(3,6)]
for lab,co in [('u',{1:1}),('u+u^3',{1:1,3:1}),('u+u^2',{1:1,2:1}),('u+u^2+2u^3',{1:1,2:1,3:2}),('u+3u^2+6u^3 (exp jet)',{1:1,2:3,3:6})]:
    print(f'{lab:24s} D6 {hilbert(D6,co,P)}  E6 {hilbert(E6,co,P)}')
