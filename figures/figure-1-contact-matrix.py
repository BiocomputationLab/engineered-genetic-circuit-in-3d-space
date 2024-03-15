from chromosome import Chromosome

import matplotlib.pyplot as pyplot
import seaborn

DATA = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data"
OUT = "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/out"

c = Chromosome(f"{DATA}/pputidakt2440.gb", f"{DATA}/6h_S46_Expression.csv", f"{DATA}/contact_wt78.abc")

f, a = c.contact_map

