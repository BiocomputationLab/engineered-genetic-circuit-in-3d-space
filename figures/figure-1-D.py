from gene3d.chromosome import Chromosome
from gene3d.insert import inserts_from_csv, filter_constructs, gradiate_colour, COLOURS
from gene3d.utils import circle_distance

from math import log, log2, exp

import matplotlib.pyplot as pyplot
from matplotlib import colormaps
from matplotlib import colors
from matplotlib import cm
import seaborn
import numpy
import pathlib

INSERT_CMAP = "Reds"

DATA = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data"
OUT = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out"
# OUT = f"{pathlib.Path(__file__).parent.parent}/out"
# DATA = f"{pathlib.Path(__file__).parent.parent}/data"


c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

construct_files = [
    "l0_stats_V1",
    "t0_stats_V0",
    "a0_stats_V1",
    "lm_stats_V1",
    "tm_stats_V0",
    "am_stats_V1"
]

constructs = inserts_from_csv([f"{DATA}/{fn}.csv" for fn in construct_files], c)
levels = {"t0": 60, "tm": 70, "a0": 100, "am": 110, "l0": 80, "lm": 90}


def add_construct(circle, construct, r, strength):
    rstart = 800 + r * construct.strand
    rfinish = 800 + (r + 10) * construct.strand
    colour = gradiate_colour(construct.colour, strength)
          
    circle.barplot(
        "chromosome",
        data = [construct.strand * 0.9],
        positions = [construct.sstart],
        width = [abs(construct.strand * 5000)],
        base_value = 0,
        raxis_range = (rstart, rfinish),
        linewidth=0.1,
        facecolor=colormaps[INSERT_CMAP](strength),
    )
    return circle

def add_rna_seq(circle, chromosome, rmin, rmax):
    v, p, w = [], [], []
    for locus, tpm in chromosome.rnaseq.items():
        try:
            gene, t, _ = chromosome[locus]
            assert tpm == t
            v.append(log2(tpm + 1) * gene.location.strand)
            p.append(gene.location.start)
            w.append(gene.location.end - gene.location.start)
        except KeyError:
            print(f"Unable to find {gene.qualifiers['locus_tag'][0]} in chromosome")
            pass

    circle.barplot("chromosome", data=v, positions=p, width=w, base_value=0, raxis_range=(rmin, rmax))
    return circle

def add_heatmap(circle, data, rmin, rmax):
    resolution = 5000
    n = len(data)
    positions = [i * resolution for i in range(n)]
    circle.heatmap(
        "chromosome",
        data=data,
        positions=positions,
        width=resolution,
        raxis_range=(rmin, rmax),
        vmin=min(data),
        vmax=max(data),
        cmap=pyplot.cm.viridis
    )
    return circle

def measurement(construct):
    baseline = next(filter(lambda x: x.name == f"{construct.name[0]}0Tn7", constructs))
    return log((construct.mean_mean / construct.mean_mean_0) / (baseline.mean_mean / baseline.mean_mean_0))

def tn7_ring():
    circle = c.circle
    tn7_hic = [log2(x + 1) for x in c.hic[c.bp_to_bin(6170441), :]]
    tn7_hic[c.bp_to_bin(6170441)] = 0
    circle = add_heatmap(circle, tn7_hic, 801, 850)
    circle.figure.savefig(f"{OUT}/tn7-ring.svg", transparent=True)
    return circle

def mean_interaction_ring():
    circle = c.circle
    open_hic = c.openness
    circle = add_heatmap(circle, open_hic, 750, 799)
    circle.figure.savefig(f"{OUT}/openness-ring.svg", transparent=True)
    return circle

def inserts_ring():
    circle = c.circle
    for code, level in levels.items():
        rstart = 800 + level
        rfinish = 800 + (level + 10)
        circle.lineplot(
            "chromosome",
            data=numpy.ones(1237),
            positions=[i * 5000 for i in range(1237)],
            raxis_range=(rstart, rfinish),
            rlim=(0, 2),
            linecolor=tuple(c / 255 for c in COLOURS[code]),
            linewidth=0.25,
        )
        rstart = 800 - level
        rfinish = 800 - (level + 10)
        circle.lineplot(
            "chromosome",
            data=numpy.ones(1237),
            positions=[i * 5000 for i in range(1237)],
            raxis_range=(rstart, rfinish),
            rlim=(0, 2),
            linecolor=tuple(c / 255 for c in COLOURS[code]),
            linewidth=0.25,
        )
    circle.figure.savefig(f"{OUT}/inserts-ring.svg", transparent=True)
    return circle

def inserts_bars():
    circle = c.circle
    minimum = min(measurement(x) for x in constructs)
    maximum = max(measurement(x) for x in constructs)
    for construct in constructs:
        strength = (measurement(construct) - minimum) / (maximum - minimum)
        print(construct.name, measurement(construct), strength)
        add_construct(circle, construct, levels[construct.name[0:2]], strength)

    circle.figure.savefig(f"{OUT}/inserts-bars.svg", transparent=True)
    return circle

def inserts_colour_bar():
    cmap = colormaps[INSERT_CMAP]
    minimum = min(measurement(x) for x in constructs)
    maximum = max(measurement(x) for x in constructs)
    norm = colors.Normalize(vmin=minimum, vmax=maximum)
    f, a = pyplot.subplots(figsize=(1, 6), layout="constrained")
    f.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=a, orientation="vertical")
    f.savefig(f"{OUT}/inserts-colour-map.svg")
    return None

def rnaseq_ring():
    circle = c.circle
    add_rna_seq(circle, c, 900, 1000)
    circle.figure.savefig(f"{OUT}/rnaseq-ring.svg", transparent=True)

def four_corners_ring():
    circle = c.circle
    bp = []
    for locus in ["PP_3333", "PP_0642", "PP_4720", "PP_2360"]:
        bp.append(c[locus][0].location.start)
    circle.scatterplot(
        "chromosome",
        data=numpy.ones(4),
        positions=bp,
        raxis_range=(749, 750),
        markersize=10,
    )
    circle.figure.savefig(f"{OUT}/four_corners_ring.svg", transparent=True)


def lib_to_tn7(libcode):
    baseline = next(filter(lambda x: x.name == f"{libcode[0]}0Tn7", constructs))
    lib = list(filter(lambda x: x.name.startswith(libcode), constructs))
    fs = [exp(measurement(x))  for x in lib]
    tn7 = [c.openness[c.bp_to_bin(x.sstart)] for x in lib]
    # tn7 = [circle_distance(x.sstart, 6170441, len(c)) for x in lib]
    return tn7, fs

def plotit():
    f, a = pyplot.subplots()
    x1, y1 = lib_to_tn7("lm")
    print(stats.linregress(x1, y1))
    x2, y2 = lib_to_tn7("tm")
    print(stats.linregress(x2, y2))
    x3, y3 = lib_to_tn7("am")
    print(stats.linregress(x3, y3))
    a.scatter(x1, y1, label="lm")
    a.scatter(x2, y2, label="tm")
    a.scatter(x3, y3, label="am")
    a.legend()
    a.set_xlabel("HiC interaction with Tn7 (higher is closer)")
    a.set_ylabel("Dynamic range versus insertion at Tn7")
    return f


if __name__ == "__main__":
    tn7_ring()
    mean_interaction_ring()
    inserts_ring()
    inserts_bars()
    inserts_colour_bar()
    four_corners_ring()
    rnaseq_ring()



    
    


