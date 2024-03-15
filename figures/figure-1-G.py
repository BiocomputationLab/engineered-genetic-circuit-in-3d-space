from gene3d.chromosome import Chromosome, openness

import pathlib
import numpy
import pandas
import matplotlib.pyplot as pyplot
import seaborn
from statistics import mean
from iced import normalization
from math import log2


DATA = f"{pathlib.Path(__file__).parent.parent}/data"
OUT = f"{pathlib.Path(__file__).parent.parent}/out"
# DATA = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data"

c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

stringdb = pandas.read_csv(f"{DATA}/stringdb_transformed.csv", ",")

def confidently_coexpressed(db, confidence_level):
    return db[(db.coexpression > confidence_level) | (db.coexpression_transferred > confidence_level)]

def distribution_hic():
    hic = c.hic
    n = hic.shape[0]
    N = n * (n + 1) // 2 - n
    distribution = numpy.zeros(N)
    i = 0
    hic = normalization.ICE_normalization(hic)
    for j in range(n-1):
        for k in range(j + 1, n):
            distribution[i] = hic[j, k]
            i = 1 + i
    return distribution

fig, axes = pyplot.subplots(1, 2)
seaborn.violinplot(
    y=confidently_coexpressed(stringdb, 700).distance_bps, ax=axes[0])
axes[0].set_ylim(0, len(c) / 2)
axes[0].set_ylabel("Distance (millions of base pairs)")

seaborn.violinplot(
    y=numpy.log2(confidently_coexpressed(stringdb, 700).distance_hic + 1), ax=axes[1])
axes[1].set_ylim(bottom=0)
axes[1].set_ylabel("Log HiC interaction strength")

fig.savefig(f"{OUT}/coexpression-violins.svg")

print(mean(confidently_coexpressed(stringdb, 700).distance_hic))
print(c.hic.mean())
