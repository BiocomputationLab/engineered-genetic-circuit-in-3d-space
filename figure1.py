from chromosome import Chromosome
from insert import inserts_from_csv
from circles import add_inserts, add_rna_seq, add_hic, add_hic_ref, get_circle, add_heatmap

import matplotlib.pyplot as pyplot
import seaborn

from pycircos import Gcircle, Garc

from math import log
from statistics import mean, stdev
from scipy import stats
import numpy

C = Chromosome(
    "/home/lewis/Data/HiC/pputidakt2440.gb",
    "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data/6h_S46_Expression.csv",
    "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data/contact_wt78.abc"
)

distance_zeros = inserts_from_csv([f"/home/lewis/Data/HiC/{x}.csv" for x in ["l0_stats_V1", "a0_stats_V1", "t0_stats_V0"]], C)
distance = inserts_from_csv([f"/home/lewis/Data/HiC/{x}.csv" for x in ["lm_stats_V1", "am_stats_V1", "tm_stats_V0"]], C)
insertions = inserts_from_csv([f"/home/lewis/Data/HiC/{x}.csv" for x in ["lm_stats_V1", "l0_stats_V1", "a0_stats_V1", "am_stats_V1", "t0_stats_V0", "tm_stats_V0"]], C)

def circle_distance(a, b, N):
    d = abs(b - a)
    if d > N / 2:
        return N - d
    else:
        return d

def openness_vs_normal_to_tn7(repressor, insertions):
    tn7 = next(filter(lambda x: x.name == f"{repressor[0]}0Tn7", insertions))
    targets = list(filter(lambda x: x.name.startswith(f"{repressor}"), insertions))
    mean_tpm = [t.mean_mean / C[t.sstart][1] for t in targets]
    openness = C.openness
    opennesses= [openness[C.bp_to_bin(t.sstart)] for t in targets] 
    return opennesses, mean_tpm

def normal_to_tn7_vs_ori_distance(repressor, insertions):
    tn7 = next(filter(lambda x: x.name == f"{repressor[0]}0Tn7", insertions))
    targets = list(filter(lambda x: x.name.startswith(f"{repressor}"), insertions))
    mean_tpm = [t.mean_mean / C[int(t.sstart)][1] for t in targets]
    distances = [circle_distance(8947, t.sstart, len(C)) for t in targets]
    return distances, mean_tpm

def plot_normal_to_tn7_vs_ori_distance(insertions):
    f, a = pyplot.subplots()
    libs = ["l0", "t0", "a0", "lm", "tm", "am"]
    xy = [normal_to_tn7_vs_ori_distance(a, insertions) for a in libs]
    for ((x, y), name) in zip(xy, libs):
        r = stats.linregress(x, y)
        a.scatter(x, y, label=name)
        a.plot([min(x), max(x)], [r.intercept + r.slope * min(x), r.intercept + r.slope * max(x)])
        print(r.intercept, r.slope, r.pvalue, min(x), max(x))

    a.set_yscale("log")
    a.set_xlabel("distance")
    a.set_ylabel("")
    a.legend()
    
    return f

def plot_openness_vs_normal_to_tn7(insertions):
    f, a = pyplot.subplots()
    libs = ["l0", "t0", "a0", "lm", "tm", "am"]
    xy = [openness_vs_normal_to_tn7(a, insertions) for a in libs]
    for ((x, y), name) in zip(xy, libs):
        r = stats.linregress(x, y)
        a.scatter(x, y, label=name)
        a.plot([min(x), max(x)], [r.intercept + r.slope * min(x), r.intercept + r.slope * max(x)])
        print(r.intercept, r.slope, r.pvalue, min(x), max(x))
        
    a.set_xlabel("distance")
    a.set_ylabel("")
    a.legend()
    
    return f

def linear_correlation_plot(x, y):
    f, a = pyplot.subplots()
    r = stats.linregress(x, y)
    print(r.intercept, r.slope, r.pvalue, min(x), max(x))
    a.set_xscale("log")
    a.scatter(x, y, s=10)
    a.plot([min(x), max(x)], [r.intercept + r.slope * min(x), r.intercept + r.slope * max(x)], c="red")
    a.set_xlabel("TPM")
    a.set_ylabel("Openness")
    return f

def linear_correlation_histogram(x, y):
    f, a = pyplot.subplots()
    r = stats.linregress(x, y)
    print(r.intercept, r.slope, r.pvalue, min(x), max(x))
    a.hist2d(numpy.log(x), y, bins=100)
    a.plot(
        [log(min(x)), log(max(x))],
        [r.intercept + r.slope * min(x), r.intercept + r.slope * max(x)],
        c="red"
    )
    return f

def family_circle(chromosome):
    circle = get_circle(chromosome)
    add_rna_seq(circle, chromosome, chromosome.forward_strand, 800, 900)
    add_rna_seq(circle, chromosome, chromosome.reverse_strand, 600, 700)
    logtn7 = numpy.log(numpy.add(chromosome.hic[1234, :], 1))
    add_heatmap(circle, chromosome.openness, 755, 795)
    add_heatmap(circle, logtn7, 705, 745)
    add_inserts(circle, distance_zeros, 590, 910, "xkcd:cyan")
    add_inserts(circle, distance, 570, 930, "xkcd:wine")
    return circle

def distance_circle(chromosome):
    circle = get_circle(chromosome)
    # add_rna_seq(circle, chromosome, chromosome.forward_strand, 800, 900)
    # add_rna_seq(circle, chromosome, chromosome.reverse_strand, 600, 700)
    logtn7 = numpy.log(numpy.add(chromosome.hic[1234, :], 1))
    # add_heatmap(circle, chromosome.openness, 755, 795)
    add_heatmap(circle, logtn7, 850, 900)
    # add_inserts(circle, distance_zeros, 590, 910, "xkcd:cyan")
    add_inserts(circle, distance, 820, 930, "xkcd:wine")
    return circle

def distance0_circle(chromosome):
    circle = get_circle(chromosome)
    # add_rna_seq(circle, chromosome, chromosome.forward_strand, 800, 900)
    # add_rna_seq(circle, chromosome, chromosome.reverse_strand, 600, 700)
    # logtn7 = numpy.log(numpy.add(chromosome.hic[1234, :], 1))
    add_heatmap(circle, chromosome.openness, 850, 900)
    # add_heatmap(circle, logtn7, 850, 900)
    add_inserts(circle, distance_zeros, 820, 930, "xkcd:cyan")
    # add_inserts(circle, distance, 820, 930, "xkcd:wine")
    return circle    

def figure1b(chromosome, filename):
    fig, ax = chromosome.contact_map
    fig.savefig(filename, dpi=900)
    return fig
    
def figure1c(chromosome, insertions, filename):
    circle = Gcircle()
    N = len(chromosome)
    circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
    circle.set_garcs()

    add_hic(circle, chromosome)
    add_rna_seq(circle, chromosome)
    add_inserts(circle, insertions)

    circle.figure.savefig(filename)
    return circle

def figure1d(chromosome, insertions, filename):
    circle = Gcircle()
    N = len(chromosome)
    circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
    circle.set_garcs()

    add_hic_ref(circle, chromosome, 1234)
    add_rna_seq(circle, chromosome)
    add_inserts(circle, insertions)
    
    circle.figure.savefig(filename)
    return circle

# figure1b(C, "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out/contact_matrix.svg")
# figure1c(C, distance_zeros, "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out/openness.svg")
# figure1d(C, distance, "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out/hic_tn7.svg")

def correlation_for_inserts(chromosome, inserts):
    transcriptions = []
    opennesses = []
    for insert in inserts:
        if insert.tpm > 0:
            transcriptions.append(insert.tpm)
            opennesses.append(insert.openness)
    return linear_correlation_plot(transcriptions, opennesses)

def correlation_for_bins(chromosome):
    transcriptions = []
    openness = chromosome.openness
    for i in range(len(openness)):
        gene, tpm, hic = chromosome[5000 * (2 * i + 1) // 2]
        transcriptions.append(tpm)
    return linear_correlation_plot(transcriptions, openness)

def zscore(x, y):
    return (x - mean(y)) / stdev(y)

def correlation_for_rna_seq(chromosome):
    transcriptions = []
    opennesses = []
    openness = chromosome.openness
    for locus in chromosome.rnaseq.keys():
        gene, tpm, hic = chromosome[locus]
        if tpm > 0.0 and abs(zscore(tpm, chromosome.rnaseq.values())) < 5:
            bp = gene.location.start if gene.location.strand == 1 else gene.location.end
            transcriptions.append(tpm)
            opennesses.append(openness[chromosome.bp_to_bin(bp)])
    return linear_correlation_plot(transcriptions, opennesses)


# def inliers(y, z):
#     return list(filter(lambda x: abs(zscore(x, y)) < z, y))

# r = stats.linregress(transcriptions, opennesses)
# f, a = pyplot.subplots()
# a.scatter(transcriptions, opennesses)
# a.plot([0, 1000], [r.intercept, r.intercept + r.slope * 1000], c="red")
# a.set_xscale("log")

