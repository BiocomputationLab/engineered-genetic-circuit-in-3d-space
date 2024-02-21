from pycircos import Garc, Gcircle
import matplotlib.pyplot as pyplot
import math
import numpy
from scipy import stats
import pandas
import random
from statistics import mean
import statistics

from data import rna_seq, contact_matrix, rna_seq_bars, libraries

sample = "gfp82"
circle = Gcircle()
rna = rna_seq(46)
contacts = contact_matrix(sample)
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
    temps = [sum(matrix[i, :]) - matrix[i, i] for i in range(matrix.shape[0])]
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=0, vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_hic_ref(circle, matrix, reference):
    temps = numpy.log2(numpy.add(matrix[reference, :], 1))
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=min(temps), vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_spots(circle, points, strand):
    rrange = (850, 860) if strand == 1 else (590, 600)
    circle.scatterplot(
        "chromosome",
        data=[0] * len(points),
        positions=points,
        markersize=20,
        rlim=(-0.1, 0.1),
        facecolor="xkcd:cyan",
        raxis_range=rrange,
    )

def add_tus(circle, df, colour):
    forward = df[(df.Direction == "forward")]
    reverse = df[(df.Direction == "reverse")]
    points = (forward.Maximum + forward.Minimum) / 2
    circle.scatterplot(
        "chromosome",
        data=[0] * len(points),
        positions=points,
        markersize=20,
        rlim=(-0.1, 0.1),
        facecolor=colour,
        raxis_range=(850, 860),
    )
    points = (reverse.Maximum + reverse.Minimum) / 2
    circle.scatterplot(
        "chromosome",
        data=[0] * len(points),
        positions=points,
        markersize=20,
        rlim=(-0.1, 0.1),
        facecolor=colour,
        raxis_range=(590, 600),
    )

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


add_hic_ref(circle, contacts, 3)
#add_hic(circle, contacts)
add_rna_seq(circle, forward, reverse)
add_library(circle, insertions)
# add_spots(circle, [(54901 + 54227) / 2, (3515484 + 3515960) / 2], 1)
# add_spots(circle, [(27599 + 27727) / 2, (2311596 + 2312855) / 2], -1)


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
    return sum(contacts[b, :]) - contacts[b, b]

def mean_openness_at(bp):
    b = int(math.floor(bp / 5000))
    return mean(contacts[b, 0:1236])

def line(x):
    return 38429.893886627586 + 1304.5248716535793 * x

def above(row):
    bp = (row.Maximum + row.Minimum) / 2
    return openness_at(bp) > line(row.TPM)

def below(row):
    bp = (row.Maximum + row.Minimum) / 2
    return openness_at(bp) < line(row.TPM)

lower = rna[[below(row) for row in rna.iloc()]]
upper = rna[[above(row) for row in rna.iloc()]]

# add_tus(circle, lower, "xkcd:wine")
# add_tus(circle, upper, "xkcd:cyan")

inserted_bins = [int(math.floor(bp / 5000)) for bp in insertions.sstart]
insertion_counts = [inserted_bins.count(i) for i in range(contacts.shape[0])]
opennesses = [sum(contacts[i, :]) - contacts[i, i] for i in range(contacts.shape[0])]
insertion_mean_openness = mean(numpy.multiply(insertion_counts, opennesses))

fig_inserts, ax_inserts = pyplot.subplots()
ax_inserts.scatter(opennesses, insertion_counts)

# fig_hist_openness, ax_hist_openness = pyplot.subplots()
# orange = range(int(max(opennesses)))
# kde_all = stats.gaussian_kde(opennesses, 0.3)
# kde_inserts = stats.gaussian_kde([openness_at(i) for i in insertions.sstart], 0.3)
# ax_hist_openness.plot(orange, kde_all(orange), lw=2, label="Whole chromosome")
# ax_hist_openness.plot(orange, kde_inserts(orange), lw=2, c="red", label="Inserts only")
# ax_hist_openness.set_xlabel("Openness")
# ax_hist_openness.set_ylabel("Density")
# ax_hist_openness.legend()

# bps = random.sample(range(-N-1, N+1), 2048)
# x = numpy.zeros(len(bps))
# for i in range(len(bps)):
#     x[i] = transcription_at(bps[i])

# x1 = numpy.log(rna.TPM + 1)
# y1 = [mean_openness_at(m) for m in (rna.Maximum + rna.Minimum) / 2]

# fig, ax = pyplot.subplots()
# ax.set_title(f"sample {sample}")
# ax.set_ylabel("Openness")
# ax.set_xlabel("TPM")
# ax.scatter(x1, y1)
# # ax.scatter([list(x1)[3338],], [list(y1)[3338],], c="red")
# r = stats.linregress(x1, y1)
# ax.plot([1, max(x1) + 1], [r.intercept, r.intercept + r.slope * math.log(max(x1) + 1)], "red")
# # ax.set_yscale('log')
# ax.set_xscale('log')

def tpm_whole_chromosome():
    x = 0
    for insert in insertions.iloc():
        L = insert.Maximum - insert.Minimum

def plot_tpm_all_versus_tpm_inserts(M):
    bps = random.sample(range(-N-1, N+1), M)
    x = numpy.zeros(len(bps))
    for i in range(len(bps)):
        x[i] = transcription_at(bps[i])
    y = numpy.zeros(len(insertions))
    for i in range(len(y)):
        bp = insertions.sstart[i]
        if insertions.sstrand[i] == "minus":
            bp = -bp
        y[i] = transcription_at(bp)

    fig, ax = pyplot.subplots()
    ax.set_ylabel("Density")
    ax.set_xlabel("TPM")
    kde_all = stats.gaussian_kde(x)
    kde_inserts = stats.gaussian_kde(y)
    orange = range(int(max(y)))
    ax.plot(orange, kde_all(orange), lw=2, label="Whole chromosome")
    ax.plot(orange, kde_inserts(orange), lw=2, c="red", label="Inserts only")
    ax.hist(x, density=True)
    ax.hist(y, density=True)
    ax.legend()
    return fig

def tpm_of_inserts(inserts):
    y = numpy.zeros(len(inserts))
    starts = list(inserts.sstart)
    strands = list(inserts.sstrand)
    for i in range(len(y)):
        bp = starts[i]
        if strands[i] == "minus":
            bp = -bp
        y[i] = math.log(transcription_at(bp) + 1)
    return y
    
def plot_tpm_inserts_histogram(inserts):
    y = numpy.zeros(len(inserts))
    starts = list(inserts.sstart)
    strands = list(inserts.sstrand)
    for i in range(len(y)):
        bp = starts[i]
        if strands[i] == "minus":
            bp = -bp
        y[i] = math.log(transcription_at(bp) + 1)

    print(mean(y))
    print(statistics.median(y))
    print(y)
    print(list(y).count(0.0))

    fig, ax = pyplot.subplots()
    ax.set_ylabel("Density")
    ax.set_xlabel("Log TPM")
    kde = stats.gaussian_kde(y)
    orange = range(int(max(y)))
    ax.plot(orange, kde(orange), lw=2)
    ax.hist(y, density=True)
    return fig, ax
    
        
d1 = pandas.read_csv("/home/lewis/sauces/python/engineered-genetic-circuits-in-3d-space/data/variants_all_positions.csv")
d1 = d1[[lib in ["5A4A", "5L14L3", "TS2dm0"] for lib in d1.Library]]

f, a = plot_tpm_inserts_histogram(d1)

covered = 0
for row in rna.iloc():
    if row.TPM > 0:
        covered = covered + row.Maximum - row.Minimum

coverage = covered / (N * 2)

