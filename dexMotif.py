##!/usr/bin/env python2.7
import os
from optparse import OptionParser
import csv
import glob
from eon import fa
from eon import dex
def check_files(taxon='mus musculus',dir = 'annotations/'):
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
	
		
def main():
	description='''
	this program tool will glast(go blasting) exons of genes given. outputting how 
	a score of how common the motifs found on the exon where present in genes within 
	go annotation categories. The assumption is that motifs important for a specific 
	gene ontology category would be present in other genes within that GO category.
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

	parser.add_option('-P','--prosite',
				action='store', type='string',
				dest='prosite',default='annotations/prosite.dat',
				help='prosite.dat file')
	
	(options, args) = parser.parse_args()
	v = options.v
	#output = options.output
	annotations = options.annotations
	prosite = options.prosite
	dexseq = options.dexFile
	if v: print 'Analyzing',dexseq
	if options.purge:
		commandOptions ='-rf'
		if v:commandOptions+='v'
		os.system('rm '+commandOptions+' annotations/*')
	
	if annotations != '':
		check_files(taxon = annotations)
	


	dexMotif = dex.dex(dexseq,prosite,v)
	dexMotif.addMotifs()

if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''
