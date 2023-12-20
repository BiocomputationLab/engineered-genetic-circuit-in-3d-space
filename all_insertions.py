from pycircos import Garc, Gcircle
import matplotlib.pyplot as pyplot
import math
import numpy
import pandas
import random

from data import rna_seq, contact_matrix, rna_seq_bars, libraries

circle = Gcircle()
rna = rna_seq(46)
contacts = contact_matrix("wt78")
insertions = pandas.read_csv("/home/lewis/sauces/python/engineered-genetic-circuits-in-3d-space/data/variants_all_positions.csv")

N = 6181873

circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
circle.set_garcs()

forward = rna[(rna.Direction == "forward")]
reverse = rna[(rna.Direction == "reverse")]

libcolours = ["xkcd:cyan", "xkcd:peach", "xkcd:dark pink", "xkcd:sky blue", "xkcd:mint green", "xkcd:wine", "xkcd:gold", "xkcd:tangerine"]

def add_rna_seq(circle, forward, reverse):
    v, p, w = rna_seq_bars(forward)
    circle.barplot(
        "chromosome",
        data=[math.log(i + 1) for i in v],
        positions=p,
        width=w,
        base_value=0.0,
        raxis_range=(760, 850),
    )
    v, p, w = rna_seq_bars(reverse)
    circle.barplot(
        "chromosome",
        data=[-math.log(i + 1) for i in v],
        positions=p,
        width=w,
        base_value=0.0,
        raxis_range=(600, 690),
    )
    return circle

def add_hic(circle, matrix):
    temps = [sum(matrix[i, :]) for i in range(matrix.shape[0])]
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_hic_ref(circle, matrix, reference):
    temps = numpy.log2(numpy.add(matrix[reference, :], 1))
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(850, 875), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_library(circle, inserts):
    forwards = inserts[(inserts.sstrand == "plus")]
    reverses = inserts[(inserts.sstrand == "minus")]
    circle.scatterplot(
        "chromosome",
        data=[0] * len(forwards),
        positions=forwards.sstart,
        markersize=20,
        rlim=(-0.1, 0.1),
        facecolor="xkcd:wine",
        raxis_range=(850, 860),
    )
    circle.scatterplot(
        "chromosome",
        data=[0] * len(reverses),
        positions=reverses.sstart,
        facecolor="xkcd:wine",
        markersize=20,
        rlim=(-0.1, 0.1),
        raxis_range=(590, 600),
    )


add_hic_ref(circle, contacts, 1235)
add_rna_seq(circle, forward, reverse)
add_library(circle, insertions)

def transcription_at(bp):
    seq_data = forward if bp > 0 else reverse
    for gene in seq_data.iloc():
        if gene.Minimum <= abs(bp) <= gene.Maximum:
            return gene.TPM
        elif gene.Minimum > abs(bp):
            break
    return 0.0

def openness_at(bp):
    b = int(math.floor(bp / 5000))
    return sum(contacts[b, :])
    

bps = random.sample(range(-N-1, N+1), 2048)
x = numpy.zeros(len(bps))
for i in range(len(bps)):
    x[i] = transcription_at(bps[i])

y = numpy.zeros(len(insertions))
for i in range(len(insertions)):
    y[i] = transcription_at(insertions.sstart[i])


