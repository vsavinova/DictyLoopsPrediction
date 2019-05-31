"""
Script to run APE threshold calculation and motifs search. 
For more information read https://github.com/VorontsovIE/sarus 

Example run:

python 02_hit_motif.py --genome dicty --transpose --background 0.388,0.112,0.112,0.388 --singleMotif "../data/motifs/autosome/pwm/jaspar/MA0139.1\ CTCF.pwm" --folderThresholds ../data/motifs_thresholds/dicty/
python 02_hit_motif.py --genome dicty --background 0.388,0.112,0.112,0.388 --folderMotifs ../data/motifs/autosome/pwm/hocomoco_11_human/ --folderThresholds ../data/motifs_thresholds/dicty/autosome/pwm/hocomoco_11_human/

python 02_hit_motif.py --genome dm3 --background 0.29,0.21,0.21,0.29 --folderMotifs ../data/motifs/pfm/ --folderThresholds ../data/motifs_thresholds/dm3/autosome/pwm/
"""

import argparse
import logging
import subprocess
import glob
import os

def call_and_check_errors(command):
	'''Function to run bash commands with subprocess module and report output.'''

	proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							shell=True, executable='/bin/bash')
	(stdout, stderr) = proc.communicate()
	logging.info("Check stdout: {}".format(stdout))
	if stderr:
		logging.info("Stderr is not empty. Might be an error in call_and_check_errors for the command: {}".format(command))
		logging.info("Check stderr: {}".format(stderr))
		return stderr
	else:
		return 0

# Command line arguments parser

parser = argparse.ArgumentParser(description='Motif search Python script.')

parser.add_argument('-g', '--genome', type=str, default='dm3', help='name of genome assembly (as in genome folder)')
parser.add_argument('--genomeFolder', type=str, default='../data/genomes', help='folder with genomes')

parser.add_argument('--apeBin', type=str, default='ape.jar', help='jre file with APE software')
parser.add_argument('--sarusBin', type=str, default='sarus.jar', help='jre file with SARUS software. Note that SARUS should be downloaded from https://github.com/VorontsovIE/sarus (support of bed output)')

parser.add_argument('-t', '--transpose', dest='transpose', action='store_true', help='Whether the pwm matrix should be transposed. False by default')
parser.set_defaults(transpose=False)

parser.add_argument('--singleMotif', type=str, default='', help='file with single motif')
parser.add_argument('--motifName', type=str, default='', help='optional motif name')
parser.add_argument('--folderMotifs', type=str, default='', help='folder with multiple motifs')
parser.add_argument('--folderThresholds', type=str, default='', help='folder with thresholds')

parser.add_argument('-T', '--thresholds', type=str, default='1e-7,0.15,1.05,mul', help='APE string with thresholds')
parser.add_argument('-i', '--ignoreThresholds', dest='ignoreThresholds', action='store_true', help='Ignore calculation of the tresholds')
parser.set_defaults(ignoreThresholds=False)
parser.add_argument('-p', '--pValue', type=float, default=1e-5, help='P-value for calculation of the score')
parser.add_argument('-b', '--background', type=str, default='0.25,0.25,0.25,0.25', help='String with nucleotides background probabilities [ACGT], comma-separated')

parser.add_argument('--outputFolder', type=str, default='tmp', help='Directory for the output. Default is "./tmp/". Note that if directory does not exist it will be created.')

parser.add_argument('--loggingLevel', type=str, default='DEBUG', help='Logging parameter: DEBUG, INFO, ERROR etc.')

args = vars(parser.parse_args()) # final dictionary with arguments

# Logging levels setup
level = logging.getLevelName(args['loggingLevel'])
logging.basicConfig(level=level)

# Calculation of thresholds for motifs hit
if not args['ignoreThresholds']:
	if len(args['singleMotif'])>0: # single motif case
		command_findThresholds = "java -cp {apeBin} ru.autosome.ape.PrecalculateThresholds {singleMotif} {folderThresholds} --background {background} --single-motif --pvalues {thresholds}".format(**args)
	else: # set of motifs case
		command_findThresholds = "java -cp {apeBin} ru.autosome.ape.PrecalculateThresholds {folderMotifs} {folderThresholds} --background {background} --pvalues {thresholds}".format(**args)
	
	if args['transpose']:
	  command_findThresholds += ' --transpose'
	
	call_and_check_errors(command_findThresholds)

# Creation of output folder
if not os.path.exists(args['outputFolder']):
	os.makedirs(args['outputFolder'])

# Running motifs search in multifasta
basic_command_runSearch = "java -cp {sarusBin} ru.autosome.SARUS {genomeFolder}/{genome}/{genome}.fa {singleMotif} {pValue:.1e} --pvalues-file {folderThresholds}/{motifName}.thr --threshold-mode pvalue --output-bed > {outputFolder}/{motifName}_{genome}_{pValue:.1e}.bed"

if len(args['singleMotif'])>0: # single motif case
	if len(args['motifName'])==0:
		args['motifName'] = '.'.join(args['singleMotif'].split('/')[-1].split('.')[:-1])
	command_runSearch = basic_command_runSearch.format(**args)
	call_and_check_errors(command_runSearch)

else: # set of motifs (all in a given folder) case
	motif_files = glob.glob(args['folderMotifs']+'/*')
	for motif_file in motif_files:
		args['singleMotif'] = motif_file
		args['motifName'] = '.'.join(args['singleMotif'].split('/')[-1].split('.')[:-1])
		command_runSearch = basic_command_runSearch.format(**args)
		call_and_check_errors(command_runSearch)
