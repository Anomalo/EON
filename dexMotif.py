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
	'''.replace('\n','').replace('\t','')
	parser = OptionParser(description=description)
	parser.add_option('-f','--file',
				action='store',type='string',
				dest='dexFile',default='',
				help='dexseq output file to analyze')

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
				dest ='v', default=False,
				help='shows you what am I thinking')
	'''
	parser.add_option('-o','--output',
				action='store', type='string',
				dest='output',default='results',
				help='directory where to save the results')
	'''
	
	(options, args) = parser.parse_args()
	v = options.v
	#output = options.output
	annotations = options.annotations
	annotaion_version = options.annotation_version
	dexseq = options.dexFile
	if v: print 'Analyzing',dexseq
	if options.purge:
		commandOptions ='-rf'
		if v:commandOptions+='v'
		os.system('rm '+commandOptions+' annotations/*')
		check_files(taxon = annotations,version =annotation_version )
	
	if annotations != '':
		check_files(taxon = annotations,version =annotation_version )
	


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
