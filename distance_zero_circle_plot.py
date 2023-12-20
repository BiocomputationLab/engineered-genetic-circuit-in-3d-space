from pycircos import Garc, Gcircle
import matplotlib.pyplot as pyplot
import math
import numpy

from data import rna_seq, contact_matrix, rna_seq_bars, libraries

lib = libraries("a0")
circle = Gcircle()
rna = rna_seq(46)
contacts = contact_matrix("wt78")

N = 6181873

circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
circle.set_garcs()

forward = rna[(rna.Direction == "forward")]
reverse = rna[(rna.Direction == "reverse")]

libcolours = ["xkcd:cyan", "xkcd:peach", "xkcd:dark pink", "xkcd:sky blue", "xkcd:mint green", "xkcd:wine", "xkcd:gold", "xkcd:tangerine"]

def add_library(circle, inserts, level):
    forwards = inserts[(inserts.sstrand == "plus")]
    reverses = inserts[(inserts.sstrand == "minus")]
    circle.scatterplot(
        "chromosome",
        data=[0] * len(forwards),
        positions=forwards.sstart,
        markersize=20,
        rlim=(-0.1, 0.1),
        facecolor=libcolours[level],
        raxis_range=(850 + 10 * level, 860 + 10 * level),
    )
    circle.scatterplot(
        "chromosome",
        data=[0] * len(reverses),
        positions=reverses.sstart,
        facecolor=libcolours[level],
        markersize=20,
        rlim=(-0.1, 0.1),
        raxis_range=(590 - 10 * level, 600 - 10 * level),
    )

def add_library_labels(circle, inserts, level):
    forwards = inserts[(inserts.sstrand == "plus")]
    reverses = inserts[(inserts.sstrand == "minus")]
    for insert in forwards.iloc():
        circle.scatterplot(
            "chromosome",
            data=[0],
            positions=[insert.sstart],
            markersize=10,
            markershape=f"${insert.Name}$",
            rlim=(-0.1, 0.1),
            facecolor=libcolours[level],
            raxis_range=(850 + 10 * level, 860 + 10 * level),
        )
    for insert in reverses.iloc():
        circle.scatterplot(
            "chromosome",
            data=[0],
            positions=[insert.sstart],
            markersize=10,
            markershape=f"${insert.Name}$",
            rlim=(-0.1, 0.1),
            facecolor=libcolours[level],
            raxis_range=(590 - 10 * level, 600 - 10 * level),
        )

def add_hic(circle, matrix):
    temps = [sum(matrix[i, :]) for i in range(matrix.shape[0])]
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

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

add_hic(circle, contacts)
add_rna_seq(circle, forward, reverse)
# add_library(circle, libraries("a0"), 0)
add_library_labels(circle, libraries("a0"), 0)
# add_library(circle, libraries("t0"), 1)
add_library_labels(circle, libraries("t0"), 1)
# add_library(circle, libraries("l0"), 2)
add_library_labels(circle, libraries("l0"), 2)
# add_library(circle, libraries("am"), 3)
# add_library(circle, libraries("tm"), 4)
# add_library(circle, libraries("ts"), 5)
# add_library(circle, libraries("lm"), 6)
    


