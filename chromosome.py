import pandas
import csv
import numpy
from Bio import SeqIO
import matplotlib.pyplot as pyplot
import seaborn
import math

def tsv_to_contact_matrix(filename):
    rows = []
    with open(filename, "r") as tsv:
            for row in csv.reader(tsv, delimiter="\t"):
                rows += [[int(row[0]), int(row[1]), float(row[3])]]

    nbins = max(row[0] for row in rows) + 1
    contacts = numpy.zeros((nbins, nbins))
    for i, j, v in rows:
        contacts[i, j] = v
        contacts[j, i] = v

    return contacts


class Chromosome:
    def __init__(self, genbank, rnaseq, hic):
        self.genbank = next(SeqIO.parse(genbank, "genbank").records)

        with open(rnaseq, "r") as rnaseq_file:
            dataframe = pandas.read_csv(rnaseq_file)
            dataframe = dataframe[(dataframe.locus_tag == dataframe.locus_tag)]
        self.rnaseq = {row.locus_tag: row.TPM for row in dataframe.iloc()}
        
        self._hic = tsv_to_contact_matrix(hic)
        self._openness = [(sum(self.hic[i, :]) - self.hic[i, i]) / (self.nbins - 1) for i in range(nbins)]

    @property
    def sequence(self):
        return self.genbank.seq

    @property
    def hic(self):
        return self._hic

    @property
    def openness(self):
        return self._openness

    @property
    def nbins(self):
        return self.hic.shape[0]

    @property
    def contact_map(self):
        f, a = pyplot.subplots()
        seaborn.heatmap(
            numpy.log(self.hic + 1),
            ax=a,
            square=True,
            cbar_kws={"label": "Log HiC interaction strength"},
            rasterized=True,
            cmap="viridis",
        )
        a.invert_yaxis()
        a.set_xticks(range(0, len(self) // 5000, 200))
        a.set_xticklabels([f"{i // 200} Mbp" for i in range(0, len(self) // 5000, 200)])
        a.set_yticks(range(0, len(self) // 5000, 200))
        a.set_yticklabels([f"{i // 200} Mbp" for i in range(0, len(self) // 5000, 200)])
        a.set_xlabel("HiC Bin")
        a.set_ylabel("HiC Bin")
        return f, a

    def intergenic(self, strand = 0):
        genic = set()
        if strand == 1 or strand == 0:
            for gene in self.forward_strand:
                genic |= set(range(gene.location.start, gene.location.end + 1))
        if strand == -1 or strand == 0:
            for gene in self.reverse_strand:
                genic |= set(range(gene.location.start, gene.location.end + 1))
        return set(range(1, len(self) + 1)) - genic

    def __len__(self):
        return len(self.sequence)

    @property
    def genes(self):
        return filter(lambda f: f.type == "gene", self.genbank.features)

    @property
    def cdss(self):
        return filter(lambda f: f.type == "CDS", self.genbank.features)

    @property
    def rrna(self):
        return filter(lambda f: f.type == "rRNA", self.genbank.features)

    @property
    def forward_strand(self):
        raw = filter(lambda f: f.type == "gene" and f.location.strand == 1, self.genbank.features)
        return sorted(raw, key=lambda f: f.location.start)

    @property
    def reverse_strand(self):
        raw = filter(lambda f: f.type == "gene" and f.location.strand == -1, self.genbank.features)
        return sorted(raw, key=lambda f: f.location.start)

    def __getitem__(self, x):
        if isinstance(x, int):
            return self.get_at_base_pair(x)
        elif isinstance(x, str):
            return self.get_from_locus(x)
        else:
            raise IndexError(f"Cannot index chromosome with a {type(x)}")

    def bp_to_bin(self, bp, res=5000):
        bp = abs(bp)
        return int(math.floor(bp / res))

    def upstream_of(self, i):
        strand = 1 if i > 0 else -1
        i = abs(i)
        genes = self.forward_strand if strand == 1 else self.reverse_strand
        j = 0
        while j < len(genes) and genes[j].location.end < i:
            j = 1 + j

        if strand == 1:
            return genes[j - 1]
        else:
            return genes[j % len(genes)]
        
    def get_at_base_pair(self, i):
        strand = 1 if i > 0 else -1
        i = abs(i)
        gene = None
        for g in self.genes:
            if i in g and g.location.strand == strand:
                gene = g

        hic_bin = self.bp_to_bin(i)
        contacts = self.hic[hic_bin, :]

        TPM = 0.0
        if gene is not None:
            TPM = self.rnaseq.get(gene.qualifiers["locus_tag"][0], 0.0)
        else:
            upstream = self.upstream_of(i * strand)
            TPM = self.rnaseq.get(upstream.qualifiers["locus_tag"][0], 0.0)
            
        return (gene, TPM, contacts)

    def get_from_locus(self, s):
        gene = None
        for g in self.genes:
            if s in g.qualifiers["locus_tag"]:
                gene = g
        TPM = self.rnaseq[s]
        contacts = self.hic[self.bp_to_bin(gene.location.start if gene.location.strand == 1 else -gene.location.end), :]

        return (gene, TPM, contacts)        
            
    
