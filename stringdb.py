from chromosome import Chromosome

import matplotlib.pyplot as pyplot

import pdb

import pandas

C = Chromosome(
    "/home/lewis/Data/HiC/pputidakt2440.gb",
    "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data/6h_S46_Expression.csv",
    "/home/lewis/src/python/engineered-genetic-circuits-in-3d-space/data/contact_wt78.abc"
)

stringdb = pandas.read_csv("/home/lewis//src/python/engineered-genetic-circuits-in-3d-space/stringdb_transformed.csv", ",")



def interaction_score(locus_a, locus_b, data):
    results = data.loc[(data["protein1"] == locus_a) & (data["protein2"] == locus_b)]
    if len(results) > 0:
        return int(results.combined_score)
    else:
        return 0

def circle_distance(a, b, N):
    d = abs(b - a)
    if d > N / 2:
        return N - d
    else:
        
        return d

def bp_distance(a, b, chromosome):
    return circle_distance(a, b, len(chromosome))

def locus_distance(locus_a, locus_b, chromosome):
    location_a = chromosome[locus_a][0].location
    location_b = chromosome[locus_b][0].location

    bp_a = location_a.start if location_a.strand == 1 else location_a.end
    bp_b = location_b.start if location_b.strand == 1 else location_b.end

    return bp_distance(int(bp_a), int(bp_b), chromosome)

def cds_distance(cds_a, cds_b, chromosome):
    location_a = cds_a.location
    location_b = cds_b.location
    bp_a = location_a.start if location_a.strand == 1 else location_a.end
    bp_b = location_b.start if location_b.strand == 1 else location_b.end

    return bp_distance(int(bp_a), int(bp_b), chromosome)

all_cds = list(C.cdss)
N = len(all_cds)
distances = []
interactions = []

def get_cds(locus, chromosome):
    for cds in chromosome.cdss:
        if locus in cds.qualifiers["locus_tag"]:
            return cds    
    
def process(row):
    locus_p = row.protein1
    locus_q = row.protein2
    print(locus_p, locus_q)
    return cds_distance(get_cds(locus_p, C), get_cds(locus_q, C), C)

def process_hic(row):
    locus_p = row.protein1
    locus_q = row.protein2
    print(locus_p, locus_q)
    p = get_cds(locus_p, C).location
    q = get_cds(locus_q, C).location
    bp_p = p.start if p.strand == 1 else p.end
    bp_q = q.start if q.strand == 1 else q.end
    bin_p = C.bp_to_bin(bp_p)
    bin_q = C.bp_to_bin(bp_q)
    return C.hic[bin_p, bin_q]
    

stringdb["distance_bps"] = stringdb.apply(process, axis=1)
stringdb["distance_hic"] = stringdb.apply(process_hic, axis=1)

# n = 0
# for row in stringdb.iloc():
#     if not n % 10000:
#         print(n)
#     distance, score = process(row)
#     distances.append(distance)
#     interactions.append(score)
#     n = 1 + n        
