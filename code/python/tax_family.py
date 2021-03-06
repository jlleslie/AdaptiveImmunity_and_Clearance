#!/usr/bin/env python
'''USAGE: 
python tax_family.py family.taxonomy formatted.family.taxonomy
'''
# This script reformats the taxonomy files output by mothur to be more readable (helps with plot labels)
# by Matt Jenior 
import sys
import re

outfile = open(sys.argv[2],'w')

with open(sys.argv[1], 'r') as tax:

	outfile.write('otu\tphylum\tfamily\tgenus\n')

	tax.readline()
	for line in tax:

		otu = line.split()[0]

		taxon = line.split()[2].split(';')
		taxon = [re.sub('\(.*?\)', '', x) for x in taxon]
		
		phylum = taxon[1]
		family = taxon[4]
		genus = taxon[5] # generate a second file for this level

		outfile.write('\t'.join([otu, phylum, family, genus + '\n']))

outfile.close()

