from pycircos import Garc, Gcircle
import matplotlib.pyplot as pyplot
import math
import numpy

from data import rna_seq, contact_matrix, rna_seq_bars, libraries

refcircle = Gcircle()

refcircle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
refcircle.set_garcs()

def add_hic(circle, matrix, reference, raxis_range):
    temps = numpy.log2(numpy.add(matrix[reference, :], 1))
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=raxis_range, vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

add_hic(refcircle, contact_matrix("wt70"), 1235, (600, 700))
add_hic(refcircle, contact_matrix("gfp74"), 1235, (700, 800))
add_hic(refcircle, contact_matrix("wt78"), 1235, (800, 900))
add_hic(refcircle, contact_matrix("gfp82"), 1235, (900, 1000))

totalcircle = Gcircle()

totalcircle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
totalcircle.set_garcs()

def add_total_hic(circle, matrix, raxis_range):
    temps = [sum(matrix[i, :]) for i in range(matrix.shape[0])]
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=raxis_range, vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

add_total_hic(totalcircle, contact_matrix("wt70"), (600, 700))
add_total_hic(totalcircle, contact_matrix("gfp74"), (700, 800))
add_total_hic(totalcircle, contact_matrix("wt78"), (800, 900))
add_total_hic(totalcircle, contact_matrix("gfp82"), (900, 1000))
