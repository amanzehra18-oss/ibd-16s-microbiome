import gzip, numpy as np, sys, pickle, time
REF=sys.argv[1]; OUT=sys.argv[2]
lut=np.full(256,-1,dtype=np.int64)
for i,c in enumerate(b'ACGT'): lut[c]=i
K=8
def kmers(seq):
    a=lut[np.frombuffer(seq,dtype=np.uint8)]
    n=len(a)-K+1
    if n<=0: return np.array([],dtype=np.int64)
    code=np.zeros(n,dtype=np.int64); bad=np.zeros(n,dtype=bool)
    for j in range(K):
        x=a[j:j+n]; bad|=(x<0); code=code*4+np.where(x<0,0,x)
    return np.unique(code[~bad])
labels={}; counts=None; glist=[]
t=time.time(); nseq=0
cur=None; seqbuf=[]
def flush():
    global nseq
    if cur is None: return
    g=labels.setdefault(cur,len(labels))
    if g>=counts.shape[0]: raise SystemExit('too many')
    km=kmers(b''.join(seqbuf).upper().replace(b'U',b'T'))
    counts[g,km]+=1; prior[km]+=1; Ng[g]+=1; nseq+=1
counts=np.zeros((5200,4**K),dtype=np.uint16); prior=np.zeros(4**K,dtype=np.int64); Ng=np.zeros(5200,dtype=np.int64)
for line in gzip.open(REF,'rb'):
    line=line.strip()
    if line.startswith(b'>'):
        flush(); cur=line[1:].decode(); seqbuf=[]
        if nseq%20000==0 and nseq: print(nseq,time.time()-t,flush=True)
    else: seqbuf.append(line)
flush()
G=len(labels); counts=counts[:G]; Ng=Ng[:G]
pickle.dump(dict(labels=labels,Ng=Ng,prior=prior,nseq=nseq),open(OUT+'.meta.pkl','wb'))
np.save(OUT+'.counts.npy',counts)
print('done',nseq,G,time.time()-t)
