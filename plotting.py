from pycircos import Garc, Gcircle
import matplotlib.pyplot as pyplot
import math
import numpy

from data import rna_seq, contact_matrix, rna_seq_bars, libraries

circle = Gcircle()
rna = rna_seq(46)
contacts = contact_matrix("wt78")

N = 6181873

circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
circle.set_garcs()

forward = rna[(rna.Direction == "forward")]
reverse = rna[(rna.Direction == "reverse")]

def add_libraries(circle, library):
    pass

def add_hic(circle, matrix, reference):
    temps = numpy.log2(numpy.add(matrix[reference, :], 1))
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(850, 875), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_overall_hic(circle, matrix):
    temps = [sum(matrix[i, :]) for i in range(matrix.shape[0])]
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(850, 875), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_rna_seq(circle, forward, reverse):
    v, p, w = rna_seq_bars(forward)
    circle.barplot(
        "chromosome",
        data=[math.log(i + 1) for i in v],
        positions=p,
        width=w,
        base_value=0.0,
        raxis_range=(910, 1000),
    )
    v, p, w = rna_seq_bars(reverse)
    circle.barplot(
        "chromosome",
        data=[-math.log(i + 1) for i in v],
        positions=p,
        width=w,
        base_value=0.0,
        raxis_range=(750, 840),
    )
    return circle


