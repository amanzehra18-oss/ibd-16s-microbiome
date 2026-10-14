import numpy as np, pickle, sys, csv
S=sys.argv[1]; FA=sys.argv[2]; OUT=sys.argv[3]
m=pickle.load(open(S+'/model.meta.pkl','rb')); counts=np.load(S+'/model.counts.npy')
labels=m['labels']; Ng=m['Ng'].astype(np.float64); nseq=m['nseq']
Pk=(m['prior']+0.5)/(nseq+1.0)
names=[None]*len(labels)
for k,v in labels.items(): names[v]=k.rstrip(';').split(';')
lut=np.full(256,-1,dtype=np.int64)
for i,c in enumerate(b'ACGT'): lut[c]=i
def kmers(seq):
    a=lut[np.frombuffer(seq,dtype=np.uint8)]; n=len(a)-7
    code=np.zeros(n,dtype=np.int64); bad=np.zeros(n,dtype=bool)
    for j in range(8):
        x=a[j:j+n]; bad|=(x<0); code=code*4+np.where(x<0,0,x)
    return np.unique(code[~bad])
comp=bytes.maketrans(b'ACGT',b'TGCA')
def rc(s): return s.translate(comp)[::-1]
seqs=[l.strip().decode().upper().encode() for l in open(FA,'rb') if not l.startswith(b'>')] if False else None
# fasta: header = sequence; sequence wrapped lines; take headers
seqs=[l[1:].strip().upper() for l in open(FA,'rb') if l.startswith(b'>')]
rng=np.random.default_rng(1)
ranks=['Kingdom','Phylum','Class','Order','Family','Genus']
anc=[[tuple(n[:r+1]) for n in names] for r in range(6)]
def scores(km):
    cnt=counts[:,km].astype(np.float64)
    return np.log((cnt+Pk[km][None,:])/(Ng[:,None]+1.0))
rows=[]; orient=[]
for i,s in enumerate(seqs):
    best=None
    for o,sq in (('F',s),('R',rc(s))):
        km=kmers(sq)
        if len(km)<10: continue
        L=scores(km); tot=L.sum(1); g=int(tot.argmax())
        if best is None or tot[g]>best[0]: best=(tot[g],g,L,o)
    if best is None: rows.append([s.decode()]+[None]*6); orient.append('NA'); continue
    sc,g,L,o=best; nk=L.shape[1]; sub=max(nk//8,1)
    idx=rng.integers(0,nk,size=(100,sub)); M=np.zeros((nk,100))
    for b in range(100): np.add.at(M[:,b],idx[b],1)
    bg=(L@M).argmax(0)
    row=[s.decode()]; ok=True
    for r in range(6):
        boot=np.mean([anc[r][x]==anc[r][g] for x in bg])*100
        if ok and boot>=50 and len(names[g])>r: row.append(names[g][r])
        else: ok=False; row.append(None)
    rows.append(row); orient.append(o)
    if i%300==0: print(i,flush=True)
with open(OUT,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['ASV']+ranks+['Orientation'])
    for r,o in zip(rows,orient): w.writerow(r+[o])
print('done')
