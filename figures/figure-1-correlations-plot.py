from chromosome import Chromosome
from insert import inserts_from_csv, filter_constructs, gradiate_colour, COLOURS
from utils import circle_distance

from math import log2, floor

import matplotlib.pyplot as pyplot
import scipy.stats as stats
import seaborn
import pdb
import numpy


DATA = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data"
OUT = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out"

c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

construct_files = ["l0_stats_V1", "t0_stats_V0", "a0_stats_V1", "lm_stats_V1", "tm_stats_V0", "am_stats_V1"]

constructs = inserts_from_csv([f"{DATA}/{fn}.csv" for fn in construct_files], c)

tpms = []
openness = []
oric_strength = []
for locus, tpm in c.rnaseq.items():
    print(f"Collecting {locus}")
    gene, t, _ = c[locus]
    assert t == tpm
    middle_of_gene = int(floor((gene.location.start + gene.location.end) / 2))
    nbin = c.bp_to_bin(middle_of_gene)
    tpms.append(tpm)
    openness.append(c.openness[nbin])
    oric_strength.append(c.hic[nbin, 1])
    print(tpm, nbin, c.openness[nbin], c.hic[nbin, 1])

f, a = pyplot.subplots()
a.scatter(tpms, openness, s=4)
a.set_xscale("log")
a.set_xlabel("Transcripts per million (TPM)")
a.set_ylabel("Mean interaction strength")
# r = stats.linregress(tpms, openness)
# a.set_yscale("log")

tpms = []
openness = []
oric_strength = []
for construct in constructs:
    print(f"Collecting {construct.name}")
    position = c[int(construct.sstart)]
    nbin = c.bp_to_bin(position.bp)
    gene = position.gene
    tpm = position.TPM

    tpms.append(tpm)
    openness.append(c.openness[nbin])
    oric_strength.append(c.hic[nbin, 1])
    print(tpm, nbin, c.openness[nbin], c.hic[nbin, 1])


