from pycircos import Garc
from math import log2

from data import rna_seq, rna_seq_bars, contact_matrix

class Insertion:
    def __init__(self, row):
        self.label = row.name
        self.length = row.length
        self.strand = 1 if row.sstrand == "plus" else -1
        self.start = row.sstart
        self.end = row.send

    def as_arc(self):
        arc = Garc(arc_id=self.label, size=self.length)
        circle.add_garc(arc)
        direction = "forward" if self.strand == 1 else "reverse"
        seq_data = rna_seq(46)
        v, p, w = rna_seq_bars(seq_data[(seq_data.Direction == direction)])
        circle.barplot(
            self.label,
            data=[math.log2(i + 1) for i in v],
            positions=p,
            width=w,
            raxis_range=(900, 1000),
        )
