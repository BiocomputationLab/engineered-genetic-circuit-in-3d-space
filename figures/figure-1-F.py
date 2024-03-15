from gene3d.chromosome import Chromosome, openness

import pathlib
import numpy
import matplotlib.pyplot as pyplot

DATA = f"{pathlib.Path(__file__).parent.parent}/data"
OUT = f"{pathlib.Path(__file__).parent.parent}/out"
# DATA = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data"

c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

x = numpy.zeros(len(c.rnaseq.keys()))
y = numpy.zeros(len(c.rnaseq.keys()))
x_specials = []
y_specials = []

i = 0
for locus, tpm in c.rnaseq.items():
    gene, t, contacts = c[locus]
    assert t == tpm
    x[i] = tpm
    bp = gene.location.start if gene.location.strand == 1 else gene.location.end
    nbin = c.bp_to_bin(bp)
    y[i] = openness(contacts, nbin)
    if locus in ["PP_3333", "PP_0642", "PP_4720", "PP_2360"]:
        x_specials.append(x[i])
        y_specials.append(y[i])
    i = 1 + i

f, a = pyplot.subplots()
a.set_xlabel("Transcripts per million (TPM)")
a.set_ylabel("Mean interaction strength")
a.scatter(x, y, s=3)
a.scatter(x_specials, y_specials, s=20)
a.set_xscale("log")
f.savefig(f"{OUT}/figure-1-F.svg")

    
print("\n".join(map(lambda x: f"{x[0], x[1]}", zip(x, y))))
