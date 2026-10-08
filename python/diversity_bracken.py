"""Alpha/beta diversity and PERMANOVA from a Bracken-derived abundance table.
Input : CSV, rows = samples named like '1Disease_SRR...' / '12Healthy_SRR...', columns = taxa.
Samples 1-11 = IBD, 12-22 = healthy (Table S1). All 22 samples must be present.
Usage : python diversity_bracken.py abundance_table.csv
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova

df = pd.read_csv(sys.argv[1], index_col=0)
num = df.index.astype(str).str.extract(r"(\d+)")[0].astype(int)
df = df.set_index(num).sort_index()                 # numeric order 1..22, never string order
assert len(df) == 22 and 20 in df.index, "expected all 22 samples (sample 20 missing?)"
group = np.where(df.index <= 11, "IBD", "Healthy")
ids = df.index.astype(str).tolist()
rel = df.div(df.sum(axis=1), axis=0)

shannon = alpha_diversity("shannon", rel.values, ids=ids, validate=False)
simpson = alpha_diversity("simpson", rel.values, ids=ids, validate=False)
alpha = pd.DataFrame({"Sample": ids, "Group": group, "Shannon": shannon.values, "Simpson": simpson.values})
alpha.to_csv("S4_alpha_diversity.csv", index=False)
for c in ("Shannon", "Simpson"):
    a, b = alpha.loc[alpha.Group == "IBD", c], alpha.loc[alpha.Group == "Healthy", c]
    print(c, "IBD median %.3f, Healthy median %.3f" % (a.median(), b.median()), mannwhitneyu(a, b, alternative="two-sided"))

bray = beta_diversity("braycurtis", rel.values, ids=ids)
pd.DataFrame(bray.data, index=ids, columns=ids).to_csv("S5_bray_curtis.csv")
coords = pcoa(bray).samples.copy(); coords["Group"] = group
coords.to_csv("S6_pcoa_coordinates.csv")
perm = permanova(bray, grouping=list(group), permutations=9999)
print(perm); perm.to_frame().T.to_csv("S7_permanova.csv", index=False)
