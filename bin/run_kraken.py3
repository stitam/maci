#!/usr/bin/python3


import argparse, os,csv

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-n",help="CPU number to use",type=int,default=50)
parser.add_argument("-g",help="Input fasta (e.g. Aci1.fna)")
parser.add_argument("-o",help="Output tsv file. Default: <GENOME_NAME>.kraken.tsv")
parser.add_argument("-r",help="Output report file. Default: <GENOME_NAME>.kraken.report.tsv")
parser.add_argument("-db_name",help="Database name",default="bacteria")
parser.add_argument("-db_path",help="Database path",default="/node8_R10/kintses_lab/databases/kraken2")
parser.add_argument("-db_version",help="Database version",default="v230502")
pargs = parser.parse_args()

if not pargs.g or not os.path.isfile(pargs.g):
	print("\33[31mPlease specify an existing fna file with -g\33[0m")
	exit()

basename = pargs.g.rsplit(".",1)[0].split("/")[-1]

if not pargs.o:
	pargs.o =  "{}.kraken.tsv".format(basename)

if not pargs.r:
	pargs.r =  "{}.kraken.report.tsv".format(basename)

	
def kraken():
	cmd = "kraken2 --db {db}_{version} --threads {threads} --output {out}.kraken.raw.tsv --report {report} {genome}".format(db=os.path.join(pargs.db_path,pargs.db_name),out=basename,report=pargs.r,genome=pargs.g,threads=pargs.n,version=pargs.db_version)
	print(cmd)
	os.system(cmd)

def smartfloat(x):
	try:
		return float(x)
	except:
		return 0
	
def get_kraken_result():
	with open(pargs.r) as f, open(pargs.o,"w") as g:
		wtr = csv.writer(g,delimiter="\t")
		header = ["assembly","kraken2_taxid","kraken2_taxon_name","kraken2_taxid_freq"]
		wtr.writerow(header)
		frequency = 0
		outrow = [basename] + ["NA"]*3
		for line in f:
			block = line.strip().split("\t")
			if block[3] != "S":
				continue
			else:
				_freq = smartfloat(block[0].strip())
				if _freq > frequency:
					outrow[1:] = [block[4].strip(),block[5].strip(),block[0].strip()]
					frequency = _freq
				elif _freq == frequency:
					for ix in [1,2]:
						outrow[ix] = "{}|{}".format(outrow[ix],block[ix+3].strip())
		wtr.writerow(outrow)

kraken()
get_kraken_result()



