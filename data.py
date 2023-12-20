import pandas
import pathlib
import csv
import numpy

def libraries(*args):
    path = pathlib.Path(__file__)
    files = [
        f"{path.parent}/data/KT2440-TS2MmeI_Library_variants.csv",
        f"{path.parent}/data/KT2440_pBLAM1-2-6A_library.csv",
        f"{path.parent}/data/KT2440-6T_Library_variants.csv",
        f"{path.parent}/data/KT2440-6L-3_Library_variants.csv",
        f"{path.parent}/data/KT2440-4T5T_Library_variants.csv",
        f"{path.parent}/data/KT2440-5A4A_Library_variants.csv",
        f"{path.parent}/data/KT5L14L3_Library_variants.csv",
        f"{path.parent}/data/KT2440-TS2dm0_Library_variants.csv",
    ]
    result = None
    for fn in files:
        with open(fn, "r") as f:
            if result is None:
                result = pandas.read_csv(f)[["Name", "length", "sstart", "send", "sstrand"]]
            else:
                df = pandas.read_csv(f)[["Name", "length", "sstart", "send", "sstrand"]]
                result = pandas.concat([result, df])
    if len(args) > 0:
        result = result[[name.startswith(args[0]) for name in result.Name]]
    return result
        
def rna_seq(sample):
    path = pathlib.Path(__file__)
    with open(f"{path.parent}/data/6h_S{sample}_Expression.csv", "r") as f:
        df = pandas.read_csv(f)
        df = df[(df.locus_tag == df.locus_tag)]
        return df

def contact_matrix(sample):
    path = pathlib.Path(__file__)
    rows = []
    with open(f"{path.parent}/data/contact_{sample}.abc", "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            rows += [[int(row[0]), int(row[1]), float(row[3])]]
    nbins = max(row[0] for row in rows)
    matrix = numpy.zeros((nbins, nbins))
    for i, j, v in rows:
        matrix[i - 1, j - 1] = v
    return matrix

def rna_seq_bars(data):
    values = []
    positions = []
    widths = []
    for row in data.iloc():
        values.append(row.TPM)
        positions.append(row.Minimum)
        widths.append(row.Length)
    return values, positions, widths

        
