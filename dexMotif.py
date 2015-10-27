##!/usr/bin/env python2.7
import os
from optparse import OptionParser
import csv
import glob
from eon import fa
from eon import dex
import sys
def err(*vars):
	sys.stderr.write(' '.join(map(str,vars))+'\n')

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
	dexMotif [options] <dexseq>

	'''.replace('\n','').replace('\t','')
	parser = OptionParser(description=description)

	parser.add_option('-A','--annotation',
				action='store', type='string',
				dest='annotations',default='',
				help='downloads and generates anotation files of defined taxon -A "mus_musculus"')

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

	parser.add_option('-B','--Bsub',
				action='store_true',
				dest ='bsub', default=False,
				help='run each dexseq in bsub')
	
	parser.add_option('-b','--Bsub_options',
				action='store', type='string',
				dest ='bsub_options', default='',
				help='options for bsubs')

	
	parser.add_option('-V','--version',
				action='store_true',
				dest ='version', default=False,
				help='prints version')
	

	(options, args) = parser.parse_args()
	v = options.v
	#output = options.output
	annotations = options.annotations
	annotation_version = options.annotation_version
	bsub = options.bsub
	bsub_options = options.bsub_options
	
	if options.version:
		err('this is version',0)
	
	if options.purge:
		commandOptions ='-rf'
		if v:commandOptions+='v'
		os.system('rm '+commandOptions+' annotations/*')
		check_files(taxon = annotations,version =annotation_version )
	
	if annotations != '':
		check_files(taxon = annotations,version =annotation_version )
	

	for dexseqArg in args:
		for dexseq in glob.glob(dexseqArg):
			if v: err('Analyzing',dexseq)
			if bsub:
				cmd = 'bsub %(bsub_options)s"python dexMotif.py %(dexseq)s"'%locals()
				if v:err( cmd)
				os.system(cmd)
			else:
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
