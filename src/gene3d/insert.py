from pandas import read_csv
from matplotlib import colormaps
from math import floor

def gradiate_colour(colour, i):
    darker_colour = [c / 10 for c in colour]
    r = darker_colour[0] + (colour[0] - darker_colour[0]) * i
    g = darker_colour[1] + (colour[1] - darker_colour[1]) * i
    b = darker_colour[2] + (colour[2] - darker_colour[2]) * i
    return [int(floor(c)) for c in [r, g, b]]


DARKER = 0.6
LIGHTER = 1
COLOURS = {
    "t0": [x * DARKER for x in [86, 180, 233]],
    "tm": [x * LIGHTER for x in [86, 180, 233]],
    "a0": [x * DARKER for x in [0, 158, 115]],
    "am": [x * LIGHTER for x in [0, 158, 115]],
    "l0": [x * DARKER for x in [204, 121, 167]],
    "lm": [x * LIGHTER for x in [204, 121, 167]],
    "TS": [x * LIGHTER for x in [204, 121, 167]],
}


def filter_constructs(constructs, code):
    return filter(lambda x: x.name.startswith(code), constructs)


class Insert:
    def __init__(self, series, chromosome):
        for x in list(series.index):
            setattr(self, x, series[x])
        self.chromosome = chromosome
        self.strand = 1 if self.sstrand =="plus" else -1
        self.sstrand = 1 if self.sstrand =="plus" else -1

    @property
    def gene(self):
        position = self.chromosome[int(self.sstart * self.sstrand)]
        return position.gene

    @property
    def tpm(self):
        position = self.chromosome[int(self.sstart * self.sstrand)]
        return position.tpm

    @property
    def openness(self):
        position = self.chromosome[int(self.sstart * self.sstrand)]
        return position.openness

    @property
    def colour(self):
        return COLOURS[self.name[0:2]]

def inserts_from_csv(fns, chromosome):
    results = []
    for fn in fns:
        insertions = read_csv(fn)
        for i in insertions.iloc():
            results.append(Insert(i, chromosome))
    return results
