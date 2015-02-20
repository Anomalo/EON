#!/usr/bin/env python2.7
import os
from optparse import OptionParser


def check_files(dir = 'annotations/',taxon='mus musculus'):
	'''
	checks, downloads and installs annotation files
	'''
	try: os.makedirs(dir)
	except OSError: pass
	#check gmt
	f = glob.glob(dir+'*.gmt')
	if f == []: 
		go.getGMT(taxon = taxon, dir = dir)
		reload(go)
	#check gtf
	f = glob.glob(dir+'*.gtf')
	if f == []: 
		import GTF
		gtf.getGTF(taxon = taxon, dir = dir)
		reload(gtf)
	#check obo
	f = glob.glob(dir+'*.obo')
	if f == []: 
		import obo
		obo.getOBO(dir = dir)
		reload(obo)
	
		
def main():
	description='''
	this program tool will glast(go blasting) exons of genes given. outputting how 
	a score of how common the motifs found on the exon where present in genes within 
	go annotation categories. The assumption is that motifs important for a specific 
	gene ontology category would be present in other genes within that GO category.
	'''.replace('\n','').replace('\t','')
	parser = OptionParser(description=description)
	parser.add_option('-G','--gene',
				action='store',type='string',
				dest='gene',default='',
				help='analyze de exons of given genes')
	parser.add_option('-A','--annotation',
				action='store', type='string',
				dest='annotations',default='',
				help='downloads and generates anotation files of defined taxon -A "mus musculus"')
	parser.add_option('-p','--purge',
				action='store_true',
				dest='purge', default=False,
				help='delete all anotation files (usefull for changing the taxon of interest)')
	
	parser.add_option('-f','--file',
				action='store',type='string',
				dest='input',default='',
				help='chose a file with gene names to analyize')
	parser.add_option('-v','--verbose',
				action='store_true',
				dest ='v', default=False,
				help='shows you what am I thinking')
	parser.add_option('-s','--wordsize',
				action='store', type='int',
				dest='ws',default=11,
				help='nmer size for blasting (default is 11)')
	parser.add_option('-o','--output',
				action='store', type='string',
				dest='output',default='results',
				help='directory where to save the results')
	parser.add_option('-B','--bsub',
				action='store_true',
				dest='bsub',default=False,
				help='sends individual gene glasting to its own bsub')
	parser.add_option('-b','--bsubOptions',
				action='store', type='string',
				dest='bsubOptions',default='',
				help='options to add to bsub in a string format')
	(options, args) = parser.parse_args()
	v = options.v
	ws=options.ws
	bsub=options.bsub
	bsuboptions=options.bsubOptions
	output = options.output
	if options.purge:
		os.system('rm -rf annotations')
	
	if options.annotations != '':
		check_files(options.annotations)
	genes = []
	if options.gene != '':
		genes =genes +  options.gene.split()
	if options.input != '':
		f = open(options.input,'r')
		genes = f.read()
		f.close()
		genes = genes + genes.split()
	
	if genes!=[]:
		if not bsub: 
			import blastgo
			g = blastgo.glast(v=v,output=output,ws=ws)
		for gene in genes:
			if bsub:
				cmd = ' '.join(['bsub',
						bsuboptions,
						'"./eon.py',
						'-G',gene,
						'-ws',str(ws),
						'-o',output,
						'"'
						])
				print cmd
				os.system(cmd)
						
						
			else:
				g.glastGene(gene)
				

if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''

