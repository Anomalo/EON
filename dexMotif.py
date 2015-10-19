##!/usr/bin/env python2.7
import os
from optparse import OptionParser
import csv
import glob
from eon import fa
from eon import dex
def check_files(taxon='mus musculus',dir = 'annotations/',version='GRCm38'):
	'''
	checks, downloads and installs annotation files
	'''
	try: os.makedirs(dir)
	except OSError: pass

	#check gtf
	f = glob.glob(dir+'*.gtf')
	if f == []: 
		from eon import gtf
		gtf.getGTF(taxon = taxon, dir = dir)
		reload(gtf)
	#
	fa.set_taxon(taxon=taxon,version=version)	
	
		
def main():
	description='''
	Takes as input a dexseq output file and enriches the loci with motifs
	overrepresented (compared with the rest of the gene)
	dexMotif [options] dexseq1 dexseq2 ...

	'''.replace('\n','').replace('\t','')
	parser = OptionParser(description=description)

	parser.add_option('-A','--annotation',
				action='store', type='string',
				dest='annotations',default='',
				help='downloads and generates anotation files of defined taxon -A "mus musculus"')

	parser.add_option('-t','--taxon_version',
				action='store', type='string',
				dest='annotation_version',default='GRCm38',
				help='defines what genome version to download, default is "GRCm38"')

	parser.add_option('-p','--purge',
				action='store_true',
				dest='purge', default=False,
				help='delete all anotation files (usefull for changing the taxon of interest)')
	
	parser.add_option('-v','--verbose',
				action='store_true',
				dest ='v', default=True,
				help='shows you what am I thinking')
	
	(options, args) = parser.parse_args()
	v = options.v
	#output = options.output
	annotations = options.annotations
	annotation_version = options.annotation_version
	if options.purge:
		commandOptions ='-rf'
		if v:commandOptions+='v'
		os.system('rm '+commandOptions+' annotations/*')
		check_files(taxon = annotations,version =annotation_version )
	
	if annotations != '':
		check_files(taxon = annotations,version =annotation_version )
	

	for dexseqArg in args:
		for dexseq in glob.glob(dexseqArg):
			if v: print 'Analyzing',dexseq
			dexMotif = dex.dex(dexseq,v,taxon=annotations, version=annotation_version)
			dexMotif.addMotifs()

if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''
