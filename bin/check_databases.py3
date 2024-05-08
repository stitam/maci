#!/usr/bin/python3


import os,argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-kraken_name",help="Kraken DB name",default="bacteria")
parser.add_argument("-kraken_version",help="Kraken DB version",default="v230502")
parser.add_argument("-kraken_path",help="Kraken DB path",default="/node8_R10/kintses_lab/databases/kraken2")
parser.add_argument("-busco_name",help="Busco DB name",default="pseudomonadales")
parser.add_argument("-busco_version",help="Busco DB version",default="odb10")
parser.add_argument("-busco_path",help="Busco DB path",default="/node8_R10/kintses_lab/databases/busco")

pargs = parser.parse_args()

def redprint(c):
	print("\33[31m\u2716 {} \33[0m".format(c))

def greenprint(c):
	print("\33[32m\u2714 {} \33[0m".format(c))

def check_databases():
	dbs = True
	##check_busco
	dirname = os.path.join(pargs.busco_path,"lineages/{}_{}".format(pargs.busco_name,pargs.busco_version))
	if os.path.isdir(dirname):
		greenprint("Busco database {} version {} found here: {}. I only checked the existance of the folder, not its content.".format(pargs.busco_name,pargs.busco_version,dirname))
	else:
		redprint("{} not found. Please download the database {} version {} from BUSCO's website.".format(dirname,pargs.busco_name,pargs.busco_version))
		dbs = False
	#check kraken
	dirname = os.path.join(pargs.kraken_path,"{}_{}".format(pargs.kraken_name,pargs.kraken_version))
	if os.path.isdir(dirname):
		greenprint("Kraken database {} version {} found here: {}. I only checked the existance of the folder, not its content.".format(pargs.kraken_name,pargs.kraken_version,dirname))
	else:
		redprint("{} not found. Please download the database {} version {}.".format(dirname,pargs.kraken_name,pargs.kraken_version))
		dbs = False
	return dbs
	


if not check_databases():
	raise Exception("Sorry, some databases were not found") 
