from math import log
from matplotlib import pyplot

def get_circle(chromosome):
    circle = Gcircle()
    N = len(chromosome)
    circle.add_garc(Garc(arc_id="chromosome", size=N, interspace=0, linewidth=0, facecolor="#ffffffff"))
    circle.set_garcs()
    return circle

def add_inserts(circle, inserts, rmin, rmax, colour):
    circle.scatterplot(
        "chromosome",
        data = [i.sstrand for i in inserts],
        positions = [i.sstart for i in inserts],
        markersize=20,
        rlim=(-1, 1),
        facecolor=colour,
        raxis_range=(rmin, rmax),
        markershape="X"
    )
    return circle

def add_rna_seq(circle, chromosome, genes, rmin, rmax):
    v, p, w = [], [], []
    for gene in genes:
        try:
            g, t, _ = chromosome[gene.qualifiers["locus_tag"][0]]
            v.append(log(t + 1) * gene.location.strand)
            p.append(g.location.start)
            w.append(g.location.end - g.location.start)
        except KeyError:
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
    
def add_hic_ref(circle, chromosome, reference):
    matrix = chromosome.hic
    temps = numpy.log2(numpy.add(matrix[reference, :], 1))
    temps[reference] = 1.0
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=min(temps), vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle

def add_hic(circle, chromosome):
    temps = chromosome.openness
    n = len(temps)
    positions = [i * 5000 for i in range(0, n)]
    circle.heatmap("chromosome", data=temps, positions=positions, width=5000, raxis_range=(700, 750), vmin=min(temps), vmax=max(temps), cmap=pyplot.cm.viridis)
    return circle
