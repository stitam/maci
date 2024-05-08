#!/usr/bin/python3

###inputs

import argparse
import csv
import os
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
	"-bvbrc_genome",
	help="BVBRC_genome.csv downlaed from https://www.bv-brc.org/"
)
parser.add_argument(
	"-bvbrc_amr",
	help="BVBRC_genome_amr.csv downloaded from https://www.bv-brc.org/"
)
parser.add_argument(
	"-assemblies",
	help="Path to file containing prediction results"
)
parser.add_argument(
	"-output",
	help="Output file path",
	default="aci_with_bvbrc.tsv"
)
pargs = parser.parse_args()

BVBRC_GENOMES = pargs.bvbrc_genome
BVBRC_PHENOTYPES = pargs.bvbrc_amr
META_ORIG = pargs.assemblies
META_OUT = pargs.output



def genome_id_list(table):
	print("Reading {}".format(table))
	genome_ids = {}
	with open(table) as f:
		rdr = csv.reader(f,delimiter=",")
		header = next(rdr)
		acc_ix = header.index("Assembly Accession")
		id_ix = header.index("Genome ID")
		for row in rdr:
			item = row[acc_ix].split(".")[0].split("_")[-1]
			if item not in genome_ids:
				genome_ids[item] = [row[id_ix]]
			else:
				genome_ids[item].append(row[id_ix])
		return genome_ids

def resistance_data(table):
	ab_list = []
	print("Reading {}".format(table))
	with open(table) as f:
		rdr = csv.reader(f,delimiter=",")
		header = next(rdr)
		id_ix = header.index("Genome ID")
		ab_ix = header.index("Antibiotic")
		phen_ix = header.index("Resistant Phenotype")
		evidence_ix = header.index("Evidence")
		resistance_data={}
		for row in rdr:
			if row[evidence_ix] != "Laboratory Method":
				continue
			if row[id_ix] not in resistance_data:
				resistance_data[row[id_ix]] = set()
			if not row[phen_ix]:
				row[phen_ix] = "Unknown"
			resistance_data[row[id_ix]].add((row[ab_ix],row[phen_ix],))
			ab_list.append(row[ab_ix])
	return [resistance_data,sorted(set(ab_list))]

def merge_resistance_data_to_metadata(meta_orig,meta_out,genome_ids,resistance,ab_list):
	print("Adding resistance data to {}".format(meta_orig))
	with open(meta_orig) as f, open(meta_out,"w") as g:
		rdr = csv.reader(f,delimiter="\t")
		wtr = csv.writer(g,delimiter="\t")
		wtr.writerow(next(rdr)+["BVBRC_ID"] + ["BVBRC:{}".format(ab) for ab in ab_list])
		for row in rdr:
			gc_base = row[0].split(".")[0].split("_")[-1]
			res_data = []
			my_genome_ids = set()
			resistance_ab = defaultdict(lambda:set())
			if gc_base in genome_ids:
				for genome_id in genome_ids[gc_base]:
					my_genome_ids.add(genome_id)
					if genome_id in resistance:
						for ab,tipus in resistance[genome_id]:
							resistance_ab[ab].add(tipus)
			wtr.writerow(row + ["|".join(sorted(my_genome_ids))] + ["|".join(sorted(resistance_ab[ab])) for ab in ab_list])
	print("Ready with {}".format(meta_out))


if __name__ == "__main__":
	genome_ids = genome_id_list(BVBRC_GENOMES)
	resistance,ab_list = resistance_data(BVBRC_PHENOTYPES)
	merge_resistance_data_to_metadata(META_ORIG,META_OUT,genome_ids,resistance,ab_list)
