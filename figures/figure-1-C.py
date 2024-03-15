import seaborn
import numpy
import pathlib

from gene3d.chromosome import Chromosome

DATA = f"{pathlib.Path(__file__).parent.parent}/data"
OUT = f"{pathlib.Path(__file__).parent.parent}/out"

c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

f, a = c.contact_map
f.savefig(f"{OUT}/figure-1-C.svg")
