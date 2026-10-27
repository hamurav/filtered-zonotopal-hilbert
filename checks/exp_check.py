# Direct exact check that the exponential cubic jet (beta/alpha^2 = 2/3) is exceptional on small simple graphs (Section 6.8).
# Usage: python3 exp_check.py ../certificates
import sys
sys.path.insert(0, sys.argv[1])
from direct_hilbert import hilbert
P=0  # exact rational arithmetic (a prime here would compute modulo it)
G={'lollipop L(3,2)':[(2,3),(3,1),(4,5),(5,1),(1,4)],
   'house':[(1,2),(2,3),(3,4),(4,1),(1,5),(5,2)],
   'kite':[(2,3),(3,1),(1,4),(4,2),(3,4),(1,5)],
   'K23':[(1,3),(3,2),(2,4),(4,1),(1,5),(5,2)]}
for name,e in G.items():
    print(name,'generic u+u^2:',hilbert(e,{1:1,2:1},P),' exp jet u+3u^2+6u^3:',hilbert(e,{1:1,2:3,3:6},P),' u+u^2+2u^3:',hilbert(e,{1:1,2:1,3:2},P),flush=True)
