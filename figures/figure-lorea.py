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

INSERT_CMAP = "Greys"

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


def add_construct(circle, construct, strength):
    rstart = 750 + 50 * construct.strand
    rfinish = 750 + 150 * construct.strand
    colour = gradiate_colour(construct.colour, strength)
          
    circle.barplot(
        "chromosome",
        data = [construct.strand],
        positions = [construct.sstart],
        width = [abs(construct.strand * 5000)],
        base_value = 0,
        raxis_range = (rstart, rfinish),
        rlim=(min(construct.strand, 0), max(construct.strand, 0)),
        linewidth=0.1,
        facecolor=colormaps[INSERT_CMAP](strength),
    )
    circle.scatterplot(
        "chromosome",
        data = [max(construct.strand, 0)],
        positions = [construct.sstart],
        raxis_range = (rstart, rfinish + construct.strand * 10),
        markershape=f"${construct.name}$",
        markersize=80,
        rlim=(min(construct.strand, 0), max(construct.strand, 0)),
        facecolor="black",
    )
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

def mean_interaction_ring():
    circle = c.circle
    open_hic = c.openness
    circle = add_heatmap(circle, open_hic, 700, 800)
    circle.figure.savefig(f"{OUT}/openness-ring-thick.svg", transparent=True)
    return circle

def toggle_ring():
    circle = c.circle
    rstart = 850
    rfinish = 860
    circle.lineplot(
        "chromosome",
        data=numpy.ones(1237),
        positions=[i * 5000 for i in range(1237)],
        raxis_range=(rstart, rfinish),
        rlim=(0, 2),
        linewidth=0.25,
    )
    rstart = 740
    rfinish = 750
    circle.lineplot(
        "chromosome",
        data=numpy.ones(1237),
        positions=[i * 5000 for i in range(1237)],
        raxis_range=(rstart, rfinish),
        rlim=(0, 2),
        linewidth=0.25,
    )
    circle.figure.savefig(f"{OUT}/inserts-ring.svg", transparent=True)
    return circle

def toggle_bars():
    circle = c.circle
    minimum = min(measurement(x) for x in constructs)
    maximum = max(measurement(x) for x in constructs)
    for construct in constructs:
        if construct.name in ["amG8a", "a0E10a", "a0C8a", "l0B5a", "lmG5a", "lmA7c", "tmE7b", "t0G4c", "tmC4e"]:
            strength = (measurement(construct) - minimum) / (maximum - minimum)
            print(construct.name, measurement(construct), strength)
            add_construct(circle, construct, strength)
    circle.figure.savefig(f"{OUT}/lorea-bars.svg", transparent=True)
    return circle

if __name__ == "__main__":
    mean_interaction_ring()
    toggle_ring()
    toggle_bars()

