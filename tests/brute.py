# independent brute-force: max clique of G_N for small N, compare with t(N)
import sys, itertools
def sqfree(n):
    d=2
    while d*d<=n:
        if n%(d*d)==0: return False
        if n%d==0: n//=d
        d+=1
    return True
N=int(sys.argv[1])
V=[n for n in range(1,N+1) if not sqfree(n*n+1)]
adj={a:set(b for b in V if b!=a and not sqfree(a*b+1)) for a in V}
best=[]
def bk(R,P,X):
    global best
    if not P and not X:
        if len(R)>len(best): best=list(R)
        return
    if len(R)+len(P)<=len(best): return
    u=max(P|X,key=lambda v:len(adj[v]&P))
    for v in list(P-adj[u]):
        bk(R|{v},P&adj[v],X&adj[v]); P=P-{v}; X=X|{v}
bk(set(),set(V),set())
t=(N-7)//25+1 if N>=7 else 0
print("N",N,"|V|",len(V),"omega",len(best),"t(N)",t, sorted(best)[:12])
