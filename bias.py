from chromosome import Chromosome
from insert import inserts_from_csv


from pandas import read_csv
from scipy import stats
import matplotlib.pyplot as pyplot


C = Chromosome(
    "/home/lewis/Data/HiC/pputidakt2440.gb",
    "/home/lewis/sauces/python/engineered-genetic-circuits-in-3d-space/data/6h_S46_Expression.csv",
    "/home/lewis/sauces/python/engineered-genetic-circuits-in-3d-space/data/contact_gfp82.abc"
)

insertions = inserts_from_csv([f"/home/lewis/Data/HiC/{x}.csv" for x in ["lm_stats_V1", "l0_stats_V1", "a0_stats_V1", "am_stats_V1", "t0_stats_V0", "tm_stats_V0"]], C)

def scatter_it(axis, libcode, x, y):
    xs = [getattr(i, x) for i in insertions if i.name.startswith(libcode)]
    ys = [getattr(i, y) for i in insertions if i.name.startswith(libcode)]
    a.scatter(xs, ys)
    return a

a0_tpm = [i.TPM for i in insertions if i.name.startswith("a0")]
l0_tpm = [i.TPM for i in insertions if i.name.startswith("l0")]
t0_tpm = [i.TPM for i in insertions if i.name.startswith("t0")]

a0_mean = [i.mean_mean / i.mean_mean_0 for i in insertions if i.name.startswith("a0")]
l0_mean = [i.mean_mean / i.mean_mean_0 for i in insertions if i.name.startswith("l0")]
t0_mean = [i.mean_mean / i.mean_mean_0 for i in insertions if i.name.startswith("t0")]

f, a = pyplot.subplots()

