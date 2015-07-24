import os
import os.path

def chr_url(chromosome, taxon="mus_musculus"):
	url = "".join(["ftp://ftp.ensembl.org/pub/release-78/fasta/", 
			taxon,
			"/dna/Mus_musculus.GRCm38.dna.chromosome.",
			chromosome,
			".fa.gz"])
	return url
def get_fa(chromosome, taxon="mus_musculus", folder="annotations/genome"):
	if not os.path.exists(folder):
		os.makedirs(folder)
	url = chr_url(chromosome, taxon)
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
	return "annotations/genome/Mus_musculus.GRCm38.dna.chromosome."+chromosome+".fa"
	
def seq_coords(chromosome, start=0, end=-1, direction='+'):
	filename = chr_filename(chromosome)
	if not os.path.isfile(filename):
		get_fa(chromosome)	
	f = open(filename, "r")
	seq = f.read().replace("\n","")
	if direction == '+':
		return seq[start:end]
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

