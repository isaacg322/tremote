import pandas as pd
import numpy as np
import sys
import gzip

global total_individuals, total_snps, snpIDs, individuals
# int total_individuals/snps
# pd.DataFrame snpIDs

def haploseparator(current_path,total_snps,**options):
	# Generates a two numpy matrices filled with 0 and 1 ('0|1' from phased haplotypes).
	# Each contain one of the phased genotype strings from each individual in the current VCF

	haplotype1_string = [] # One for each phased genotype string
	haplotype2_string = []

	for indID in individuals:
		for i in range(0,total_snps):
			haplotype1.append(dataframed_vcf[indID][i][0])
			# A list of the first strand of haplotype for all the individuals
			# in the dataset
	for indID in individuals:
		for i in range(0,total_snps):
			haplotype2.append(dataframed_vcf[indID][i][2])

	h1_str_list = [haplotype1_string[x:x+total_snps] for x in range(0, len(haplotype1_string),total_snps)]
	h1_final = [] # Divide each original list by the number of snps in the dataset, creating
	# a new list for each individual

	h2_str_list = [haplotype2_string[x:x+total_snps] for x in range(0, len(haplotype2_string),total_snps)]
	h2_final = []

	for list_i in h1_str_list:
		list_i = map(int, list_i) # Convert strings to integers to obtain a np.matrix
		h1_final.append(list_i) # Each row in the 'final' lists is an individual
	for list_i in h2_str_list:
		list_i = map(int, list_i)
		h2_final.append(list_i)

	h1_matrix = np.matrix(h1_final)
	h2_matrix = np.matrix(h2_final)

	if options.get('haplotype') == 'h1':
		return(h1_matrix)
	if options.get('haplotype') == 'h2':
		return(h2_matrix)

	# The resulting matrices are used by index_annot() and almost_haplotypes() functions

def dict_generator(snpIDs,total_snps):
	# Generate a dictionary based on ALT nucleotides and position indexes ({POS1:A; POS2:G})

	values = (vcf_dataframe['ALT']).tolist()
	keys = range(0,total_snps)
	dictionary = dict(zip(keys,values))
	dictionary['-'] = 'N' # Adds the value and key to be used when a reference allele is found
	return(dictionary)

def index_annot(haplotype_matrix,total_individuals,dictionary):
	# Annotates if each matrix sites as its corresponding nt as annotated in the current dictionary
	# This and the previous function require you to select ONLY ONE ALTERNATIVE ALLELE!
	# There's a normalization of sites to make them bi-allelic

	temporal_indexes = [[]for x in xrange(0,total_individuals)]

	for dlist_index,sublists in enumerate(haplotype_matrix.tolist(), start= 0):
		for element_index,i in enumerate(sublists, start=0):
			if i == 1:
				temporal_indexes[dlist_index].append(element_index)# A way to
				# store indexes, counts every element inside a list and annotate
				# the count when it finds a 1.
	# Above iteration returns a list of lists with indexes were a match was found
			elif i == 0:
				temporal_indexes[dlist_index].append('-')

	translated_indexes = [[]for x in xrange(0,total_individuals)]

	for list_index, lists1 in enumerate(temporal_indexes, start= 0):
		if not lists1:
			pass
		else:
			for i in lists1:
				current_element = chr_dictionary[i]
				translated_indexes[list_index].append(current_element)
	# Above iteration returns a translated version  of the rsID in the list of
	# lists generated at the previous iteration.

	return(translated_indexes)

	# Next function is made to manage the presence of '-' in dictionaries
	# Returns a tuple containing np.arrays with unique string haplotype
	# configurations and counts for its ocurrences

def almost_haplotypes(haplotype_matrix,total_individuals,dictionary):

	annotate = index_annot(haplotype_matrix,total_individuals,dictionary)

	temp_list = [[] for x in xrange(0,total_individuals)]

	for num_list,ind_list in enumerate(annotate, start = 0):
		temp_list[num_list].append(''.join(ind_list))

	raw_report = np.unique(np.array(temp_list), return_counts = True)
	most_common_haplotypes = raw_report[0].tolist()
	most_common_counts = raw_report[1].tolist()
	common_haplotypes = sorted(zip(most_common_haplotypes,most_common_counts),
	key=lambda seq_string: seq_string[1],reverse=True)[:3] #Returns counts for ONLY the three most common haplotypes!

	return(common_haplotypes)

def each_chromosome(current_path):

	# VCF file parsing
	with gzip.open(current_path) as f:
		for line in f:
			line = line.partition('#')[0] # Remove vcf '#INFO' lines
			line = line.rstrip()
			if line != '':
				dataframed_vcf.append(line)
	dataframed_vcf = [x.split("\t") for x in dataframed_vcf]
	dataframed_vcf = pd.DataFrame(dataframed_vcf) # Working file
	print "VCF oppened as dataframe succesfully!"
	snpIDs = dataframed_vcf.iloc[:,5:6]
	snpIDs.columns = ['ALT'] # To dict_generator() EXPECTS TO HAVE A SINGLE NT! --> Cambia el parseo de esto en dict_generator()
	# para que no se muera cuando encuentra varias letras, haz que eso lo interprete como un '1' osea un sitio alternativo.

	dataframed_vcf = dataframed_vcf.iloc[:,9:] # Show only those columns with haplotype info
	individuals = dataframed_vcf.columns.values # Each individual ID in dataframed_vcf

	#Generic info for the VCF in the current iteration
	total_snps = len(dataframed_vcf)
	total_individuals = len(individuals)

	current_dict = dict_generator(snpIDs,total_snps)
	print "Haplotype dictionary generated succesfully!"

	haplo1 = haploseparator(dataframed_vcf,total_snps,
	haplotype = 'h1') # Return first phased string from haploseparator() as a np.matrix
	haplo2 = haploseparator(dataframed_vcf,total_snps,
	haplotype = 'h2') # Return second phased string from haploseparator() as a np.matrix

	translated_haplo1 = almost_haplotypes(haplo1,total_individuals,current_dict)
	translated_haplo2 = almost_haplotypes(haplo2,total_individuals,current_dict)

## Yo mejoraría desde aquí, que el output llene una tabla de cuentas de halotipo
## por región y te la guarde en un archivo, pensando en que probablemente tendrás como
## mil VCFs
	print """
	*These are the most common haplotype configurations in the
	first string of haplotype for the """,current_path,"""*
	""", translated_haplo1
	print """
	*These are the most common haplotype configurations in the
	second string of haplotype for the """,current_path,""" region*
	""", translated_haplo2



print """This program disects and counts haplotypes from phased genotype chains in
VCF files"""

start_path = raw_input("""Enter the directory containing the gzipped VCF's
to be analyzed: """)
print "Any unzipped VCF file will be ignored..."

paths = [] # Store VCF files in a given directory

for vcf_file in os.listdir(start_path):
	if vcf_file[-7:] != ".vcf.gz":
		continue
	else:
		path = os.path.join(start_path, vcf_file)
		paths.append(path)

print len(paths), " files were indexed succesfully..."
checkpoint_one = raw_input("Continue? (Y/N) ")
if checkpoint_one.lower() == "y" | checkpoint_one == "":
	print "Initializing VCF processing...\nThis may take a while..."
	for current_path in paths:
		print "Working on ", current_path
		each_chromosome(current_path)
else:
	print "*Terminating program*"
