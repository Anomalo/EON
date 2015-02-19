#!/usr/bin/env python2.7

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
	#check gob
	f = glob.glob(dir+'*.gob')
	if f == []: 
		blastgo.genGOB()

	
		
def main():
	usage = 'usage: %prog [options] arg1 arg2'
	parser = OptionParser()
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
				help='delete all anotation files (usefull for changing the taxon of interest')
	
	parser.add_option('-f','--file',
				action='store',type='string',
				dest='input',default='',
				help='chose a file with gene names to analyize')
	parser.add_option('-v','--verbose',
				action='store_true',
				dest ='v', default=False,
				help='shows you what am I thinking')
	parser.add_option('-o','--output',
				action='store', type='string',
				dest='output',default='results',
				help='directory where to save the results')
				
	(options, args) = parser.parse_args()
	v = options.v
	output = options.output
	if options.purge: os.system('rm -rf annotations')
	
	modules = ['glob','go','blastgo','fa','os']
	for module in modules:
		print 'importing',module
		exec('import '+module)
		
	if options.annotations != '':
		check_files(options.annotations)
	
	if options.gene != '':
		genes = options.gene.split()
		g = blastgo.glast()
		for gene in genes:
			g.glastGene(gene, v=v,output=output) 
	if option.input != '':
		f = open(option.input,'r')
		genes = f.read()
		f.close()
		genes = genes.split()
		g = blastgo.glast()
		for gene in genes:
			g.glastGene(gene,v=v,output=output)

if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''

