#!/bin/bash

# hp9QEOt8wQfb8NjFEMxS8yyA7Cc=

#SBATCH --job-name=HICKTWT
#SBATCH --output=./%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --partition=fast
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lewis.grozinger@upm.es

REFGENOME=pputida.fa
READS=($1 $2)
ENZYME="Sau3AI MluCI"
RES=5000
FILTERS=(1 2 3 9)
WORKDIR=hicresults-from-${READS[0]}-${READS[1]}
[! -d "${WORKDIR}" ] && mkdir -p ${WORKDIR};

MAPPER=/gem3-mapper/bin/gem-mapper
CPUS=4

function jobids {
    singularity exec hic-pipe tadbit describe --workdir ${WORKDIR} -t 2 -s Id | grep -Eo '[0-9]+'
}

function regular_mapping {
    echo "Mapping of read $1"
    echo "Using reference genome ${REFGENOME} indexed at ${REFGENOME}.gem"
    echo "Using read at ${READS[$1 - 1]}"
    echo "With restriction enzyme(s) ${ENZYME}"
    singularity exec hic-pipe \
	    tadbit map --workdir ${WORKDIR} --genome ${REFGENOME} --index ${REFGENOME}.gem \
	    --fastq ${READS[$1 - 1]} --read $1 --renz ${ENZYME} \
	    --cpus ${CPUS} --mapper_binary ${MAPPER}
}

function gem_indexing {
    echo "Indexing of reference genome ${REFGENOME}"
    singularity exec hic-pipe /gem3-mapper/bin/gem-indexer --input ${REFGENOME} --output ${REFGENOME}
}

function parse_maps {
    echo "from pytadbit.parsers.genome_parser import parse_fasta; parse_fasta(\"${REFGENOME}\")" | singularity exec hic-pipe python3
    singularity exec hic-pipe tadbit parse --workdir ${WORKDIR} --jobids 1 2 --genome ${REFGENOME} --compress_input
}

function filter_reads {
    echo "Filtering reads"
    echo "Applying filters ${FILTERS[*]}"
    singularity exec hic-pipe tadbit filter --workdir ${WORKDIR} --cpus ${CPUS} --apply ${FILTERS[*]} --compress_input --valid --clean
}

function normalise {
    echo "Normalising data"
    singularity exec hic-pipe tadbit normalize --workdir ${WORKDIR} --resolution ${RES} --cpus ${CPUS} --min_count 250
}

function bin_to_matrix {
    echo "Binning normalised data into contact matrix"
    singularity exec hic-pipe tadbit bin --workdir ${WORKDIR} --resolution ${RES} --cpus ${CPUS} --plot --norm raw norm
}

function find_tads {
    echo "Finding TADs"
    singularity exec hic-pipe tadbit segment --workdir ${WORKDIR} --resolution ${RES} --only_tads --cpu ${CPUS}
}

function find_compartments {
    echo "Finding compartments"
    singularity exec hic-pipe tadbit segment --workdir ${WORKDIR} --resolution ${RES} --only_compartments --cpu ${CPUS}
}

# BEGIN

# ml load Singularity

gem_indexing; 

regular_mapping 1; 

regular_mapping 2; 

parse_maps; 

filter_reads;

normalise;

bin_to_matrix;

find_tads;

# END
