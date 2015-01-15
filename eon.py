modules = ['glob','gtf','go','blastgo','fa','blast']
for module in modules:
	print 'importing',module
	exec('import '+module)

def eonExon(exon, taxon=''):
	chromosome, start, end = (gtf.getTranscriptCoords(exon))
	sequence = fa.seq_coords(chromosome, start, end)
	
	print 'blasting',exon
	#blast		
	b = blast.blast(sequence, program = 'blastx', filter = taxon)
	return b

def eonGene(gene, taxon=''):	
	exon_names = gtf.transcriptNames(gene)
	exons = []
	out = {}
	for exon in exon_names:
		out.update({exon:eonExon(exon, taxon)})
	return out

def check_files():
	#check gmt
	files = glob.glob('annotations/*.gmt')
	if len(files)==0:
		print 'missing gmt file in annotations/'
		return None
	else: 
		gmt = files[0]
	#check gob
	files = glob.glob('annotations/*.gob')
	if len(files)==0:
		print 'missing gob file in annotations/ \n genetating from',gmt
		blastgo.genGOB(gmt)
		return None
	gob = glob.glob('annotations/*.gob')[0]
	#check gib	
	files = glob.glob('annotations/*.gib')
	if len(files)==0:
		print 'missing gib file in annotations/ \n genetating from', gob
		blastgo.genGIB(gob)
		return None
	

if __name__ == '__main__': 
	#print eonExon('Xkr4-001')
	check_files()
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e



