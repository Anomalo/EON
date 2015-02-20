import os

def getOBO(dir = 'annotations/'):
	url ='http://geneontology.org/ontology/go.obo'
	os.system(' '.join([	'wget',
				url,
				'-P',
				dir]))
getOBO()
