import os
import os.path
from glob import glob
import sys
_taxon="Mus_musculus"
_version='GRCm38'
_annotationDir='annotations'
def set_taxon(taxon="Mus_musculus",version='GRCm38',annotationDir='annotations'):
	global _taxon
	global _version
	global _annotationDir
	_taxon = taxon
	_version = version
	_annotationDir = annotationDir

def chr_url(chromosome):
	taxon = _taxon
	version = _version
	taxon_lower = taxon.lower()
	url ="ftp://ftp.ensembl.org/pub/release-78/fasta/%(taxon_lower)s/dna/%(taxon)s.%(version)s.dna.chromosome.%(chromosome)s.fa.gz"%locals()
	return url

def get_fa(chromosome, folder="annotations/genome"):
	if not os.path.exists(folder):
		os.makedirs(folder)
	url = chr_url(chromosome)
	cmd = " ".join(["wget",
			url,
			"-P",
			folder])
	os.system(cmd)
	filename = url.split("/")[-1]
	cmd = " ".join(["gunzip",
			folder+"/"+filename])
	os.system(cmd)
	'''
	cmd = " ".join(["rm",
			folder+"/"+filename])
	os.system(cmd)
	'''

def chr_filename(chromosome):
	'''
	returns the filename for the chromosome genome
	'''
	annotationDir =  _annotationDir
	taxon=_taxon
	version = _version
	file = (("%(annotationDir)s/genome/%(taxon)s.%(version)s.dna.chromosome.%(chromosome)s.fa")%locals()).replace('//','/')
	return file

def seq_coords(chromosome, start=0, end=-1, direction='+'):
	'''
	returns the genomic sequence
	'''
	filename = chr_filename(chromosome)
	#if filename == '':
	#	get_fa(chromosome,folder = _annotationDir+'/genome')
	#filename = chr_filename(chromosome)
	try:
		f = open(filename, "r")
	except :
		get_fa(chromosome,folder = _annotationDir+'/genome')
		try:
			f = open(filename, "r")
		except:
			raise("failed to download or find the genome for taxon "+_taxon+" version "+_version+".\n please place the individual dna for each chromosome in .gz format in the a folder namded 'genome' in the annotation folder")

	seq = f.read().replace("\n","")
	start,end = int(start),int(end)
	if direction == '+':
		seq = seq[start:end]
	else:
		seqReversed = seq[start:end][::-1]
		convert = {	'A':'T',
				'C':'G',
				'G':'C',
				'T':'A',
				'N':'N'}
		seq = ''
		for i in seqReversed:
			seq+=convert[i]
	return seq

def seqs_coords(coordinates):
	'''
	return a list of sequences based on a list of coords tuples [(chromosome, start,
	end, strand),...]
	'''
	if coordinates==None:return 'N'*10
	seqs = []
	for chromosome, start, end, strand in coordinates:
		seqs.append(seq_coords(chromosome, start,end, strand))
	return seqs
