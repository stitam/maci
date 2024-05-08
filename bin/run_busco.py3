#!/usr/bin/python3


import argparse, os,csv

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-n",help="CPU number to use",type=int,default=50)
parser.add_argument("-g",help="Input fasta (e.g. Aci1.fna)")
parser.add_argument("-o",help="Output tsv file. Default: <pargs.g>.busco.tsv")
parser.add_argument("-db_name",help="Database name",default="pseudomonadales")
parser.add_argument("-db_path",help="Database path",default="/node8_R10/kintses_lab/databases/busco")
parser.add_argument("-db_version",help="Database version",default="odb10")
pargs = parser.parse_args()

if not pargs.g or not os.path.isfile(pargs.g):
	print("\33[31mPlease specify an existing fna file with -g\33[0m")
	exit()

basename = pargs.g.rsplit(".",1)[0].split("/")[-1]

if not pargs.o:
	pargs.o =  "{}.busco.tsv".format(basename)

def busco():
	os.system("export NUMEXPR_MAX_THREADS=8 && busco -m geno -l {lineage} -f -o {out} -q -i {genome} -c {threads} --datasets_version {db_version} --download_path {db_path} --offline".format(lineage=pargs.db_name,out=basename,genome=pargs.g,threads=pargs.n,db_version=pargs.db_version,db_path=pargs.db_path))

def convert_busco_to_tsv():
	with open("{base}/short_summary.specific.{db}_{version}.{base}.txt".format(base=basename,db=pargs.db_name,version=pargs.db_version)) as f, open(pargs.o,"w") as g:
		wtr = csv.writer(g,delimiter="\t")
		header = ["busco_complete","busco_single_copy","busco_multi_copy","busco_fragmented","busco_missing"]
		real_header = ["assembly"]
		for item in header:
			real_header.append(item)
			real_header.append("{}_perc".format(item))
		wtr.writerow(real_header)
		row = [basename]
		numbers = []
		percentages = []
		read = False
		for line in f:
			line = line.strip()
			if not line:
				continue
			if line.startswith("***"):
				read = True
				continue
			if line.startswith("Assembly"):
				read = False
			if read:
				if line[0] in "0123456789":
					numbers.append(line.split()[0])
				else:
					total = line.split(":")[-1]
					block = line.split(":")
					for item in block[1:-1]:
						percentages.append(item.split("%")[0])
		for i,j in zip(numbers,percentages):
			row.append(i)
			row.append(j)
		wtr.writerow(row)


busco()
convert_busco_to_tsv()
