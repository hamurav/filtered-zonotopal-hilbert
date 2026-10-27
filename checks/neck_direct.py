# Example 4.17: direct values of H_n for the necks theta*_n theta, K13*_n K13 (joined at the
# centers) and K3*_n K3 at f = u+u^2, and L(H_n) of (4.18), to compare with the neck constants
# of Theorem 4.16. Ranks modulo the prime 2^31-1, so this is a consistency check only.
# Usage: python3 neck_direct.py theta|K13|K3 NMAX   (run from checks/)
import sys, time
sys.path.insert(0,'../certificates')
from direct_hilbert import hilbert
def neck(e1,e2,n):
    off=max(x for e in e1 for x in e)+1
    e2s=[(a+off,b+off) for a,b in e2]
    # path from vertex 0 (of G') to vertex off (of G'') of length n
    path=[]; prev=0; nxt=2*off+10
    for i in range(n-1):
        path.append((prev,nxt)); prev=nxt; nxt+=1
    path.append((prev,off))
    return e1+e2s+path
theta=[(0,1),(0,1),(0,1)]
K13=[(0,1),(0,2),(0,3)]
K3=[(0,1),(1,2),(0,2)]
which={'theta':theta,'K13':K13,'K3':K3}[sys.argv[1]]
H={}
for n in range(2,int(sys.argv[2])+1):
    t=time.monotonic(); H[n]=hilbert(neck(which,which,n),{1:1,2:1},p=2147483647); print(n,H[n],f'{time.monotonic()-t:.0f}s',flush=True)
def L(n):
    m=max(len(H[k]) for k in range(n-3,n+1))+2
    g=lambda k,i: H[k][i] if 0<=i<len(H[k]) else 0
    return [g(n,i)-g(n-1,i)-3*g(n-2,i-1)+2*g(n-3,i-1) for i in range(m)]
for n in range(5,max(H)+1): print('L(H_%d)='%n, L(n))
