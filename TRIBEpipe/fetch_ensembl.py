#!/usr/bin/env python
"""download ensembl files for TRIBEpipe"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"


import os, sys, logging
import urllib.request
import pathlib

dd = {'BDGP6': {
		'gtf': 'ftp://ftp.ensembl.org/pub/release-92/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.92.gtf.gz',
		'dna': 'ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz',
	   },
	  'GRCh38': {
	   	'gtf': 'ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz',
	   	'dna': 'ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz',
	   },
	  'GRCm38': {
	  	'gtf': 'ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz',
	  	'dna': 'ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz',
	  },
	  'BDGP5': {
	  	'gtf': 'ftp://ftp.ensembl.org/pub/release-75/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.75.gtf.gz',
	  	'dna': 'ftp://ftp.ensembl.org/pub/release-75/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.75.dna_sm.toplevel.fa.gz',
	  },
	  'GRCh37': {
	  	'gtf': 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz',
	  	'dna': 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz'
	  }
}


logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)


def url_is_ok(url):
	"""Check the url is ok, accesable"""
	r = urllib.request.urlopen(url).getcode()
	if r == 200 or r is None: # https, ftp
		return True
	else:
		return False



def download_file(url, save_path = None, overwrite = False):
	"""download and save file"""
	f_name = os.path.basename(url)
	save_path = './' if save_path is None else save_path
	f_file = os.path.join(save_path, f_name)
	if os.path.exists(f_file) and overwrite is False:
		logging.info('file exists, skipping download:' + f_file)
	else:
		if url_is_ok(url):
			logging.info('start downloading file:' + f_file)
			with open(f_file, 'wb') as fo:
				fo.write(urllib.request.urlopen(url).read())
			logging.info('ok, download file:' + f_file)
		else:
			logging.error('bad url: ' + url)



def fetch_url(species = 'dm3'):
	"""ftp"""

	if species in ['dm3', 'BDGP5','fruitfly']:
		sp = 'BDGP5'
	elif species in ['dm6', 'BDGP6']:
		sp = 'BDGP6'
	elif species in ['hg38', 'GRCh38', 'human']:
		sp = 'GRCh38'
	elif species in ['hg19', 'GRCh37']:
		sp = 'GRCh37'
	elif species in ['mm10', 'GRCm38', 'mouse']:
		sp = 'GRCm38'
	else:
		logging.error('unknown species:' + species)
		return None
	## url
	url_gtf = dd[sp]['gtf']
	url_dna = dd[sp]['dna']
	return [url_gtf, url_dna]


def main():
	if(len(sys.argv) < 3):
		sys.exit('Usage: fetch_ensembl.py <dm3> <gtf> <out_path>')
	species, group = sys.argv[1:3]
	out_path = str(pathlib.Path.home()) if sys.argv[3] is '' else sys.argv[3]
	out_path = os.path.join(out_path, 'data', 'genome', species)
	url_gtf, url_dna = fetch_url(species)
	if group.lower() == 'all' or group.lower() == 'gtf':
		path_gtf = os.path.join(out_path, 'annotation_and_repeats')
		if not os.path.exists(path_gtf):
			os.makedirs(path_gtf)
		download_file(url_gtf, save_path = path_gtf)
	elif group.lower() == 'all' or group.lower() == 'dna':
		path_dna = os.path.join(out_path, 'bigZips')
		if not os.path.exists(path_dna):
			os.makedirs(path_dna)
		download_file(url_dna, save_path = path_dna)
	else:
		logging.error('unknown type <all|gtf|dna>:' + group)
	
if __name__ == '__main__':
	main()


## EOF