
import pandas as pd
import numpy as np
import sys

# This script is made to be a test version on the 532 individuals from the
# 1000 Genomes Project.
# Running this test for bigger datasets is not recomended.



chr3st_vcf = pd.read_table('~/Escritorio/interest_snps_chr3st.recode.vcf',
sep='\t')
chr3nd_vcf = pd.read_table('~/Escritorio/interest_snps_chr3nd.recode.vcf',
sep='\t')
chr5_vcf = pd.read_table('~/Escritorio/interest_snps_chr5.recode.vcf',
sep='\t')
chr7_vcf = pd.read_table('~/Escritorio/interest_snps_chr7.recode.vcf',
sep='\t')
chr20_vcf = pd.read_table('~/Escritorio/interest_snps_chr20.recode.vcf',
sep='\t')

# Full VCF for all regions in the case individuals from the original stage 2
# meta-analysis
all_cases = pd.read_table('~/Escritorio/interest_snps_chr3st.recode.vcf',
sep='\t')

def haploseparator(vcf_dataframe,number_of_snps,**options):

    nsnps = number_of_snps
    df = vcf_dataframe
    subs = df.iloc[:,9:] # Show only those columns and rows with haplotype info
    individuals = subs.columns.values

    haplotype1 = []

    haplotype2 = []

    # A list of the first strand of haplotype for all the individuals in the
    # dataset
    for column in individuals:
        for i in range(0,nsnps):
            haplotype1.append(subs[column][i][0])
            # A list of the second strand of haplotype for all the individuals
            # in the dataset
    for column in individuals:
        for i in range(0,nsnps):
            haplotype2.append(subs[column][i][2])
            # Split the list in list by the number of snps of interest

    h1_str_list = [haplotype1[x:x+nsnps] for x in range(0, len(haplotype1),nsnps)]
    h1_final = [] # Each list is an individual

    h2_str_list = [haplotype2[x:x+nsnps] for x in range(0, len(haplotype2),nsnps)]
    h2_final = [] # Each list is an individual

    for list_i in h1_str_list:
        list_i = map(int, list_i)
        h1_final.append(list_i)
    for list_i in h2_str_list:
        list_i = map(int, list_i)
        h2_final.append(list_i)

    h1_matrix = np.matrix(h1_final)
    h2_matrix = np.matrix(h2_final)

    if options.get('haplotype') == 'h1':
        return h1_matrix
    if options.get('haplotype') == 'h2':
        return h2_matrix


# Matrices for each interest region
# chr3st
h1chr3st = haploseparator(chr3st_vcf,len(chr3st_vcf),haplotype = 'h1')
h2chr3st = haploseparator(chr3st_vcf,len(chr3st_vcf),haplotype = 'h2')

# chr3nd
h1chr3nd = haploseparator(chr3nd_vcf,len(chr3nd_vcf),haplotype = 'h1')
h2chr3nd = haploseparator(chr3nd_vcf,len(chr3nd_vcf),haplotype = 'h2')

# chr5
h1chr5 = haploseparator(chr5_vcf,len(chr5_vcf),haplotype = 'h1')
h2chr5 = haploseparator(chr5_vcf,len(chr5_vcf),haplotype = 'h2')

# chr7
h1chr7 = haploseparator(chr7_vcf,len(chr7_vcf),haplotype = 'h1')
h2chr7 = haploseparator(chr7_vcf,len(chr7_vcf),haplotype = 'h2')

# chr20
h1chr20 = haploseparator(chr20_vcf,len(chr20_vcf),haplotype = 'h1')
h2chr20 = haploseparator(chr20_vcf,len(chr20_vcf),haplotype = 'h2')

# Full VCF
haplotype_1 = haploseparator(all_cases,len(all_cases),haplotype = 'h1')
haplotype_2 = haploseparator(all_cases,len(all_cases),haplotype = 'h2')

def dict_generator(vcf_dataframe,max_snps,**options):
    # Generate a dictionary based on rsIDS and position indexes
    if options.get('value') ==  'rsID':
        values = (vcf_dataframe['ID']).tolist()
    # Generate a dictionary based on alternative nucleotides and
    # position indexes
    if options.get('value') == 'nucleotide':
        values = (vcf_dataframe['ALT']).tolist()

    keys = range(0,max_snps)
    dictionary = dict(zip(keys,values))
    dictionary['-'] = 'N'
    return(dictionary)

# Nucleotide based dictionary generation for all VCF's
chr3st_dict= dict_generator(chr3st_vcf,len(chr3st_vcf), value = 'nucleotide')
chr3nd_dict= dict_generator(chr3nd_vcf,len(chr3nd_vcf), value = 'nucleotide')
chr5_dict= dict_generator(chr5_vcf,len(chr5_vcf), value = 'nucleotide')
chr7_dict= dict_generator(chr7_vcf,len(chr7_vcf), value = 'nucleotide')
chr20_dict= dict_generator(chr20_vcf,len(chr20_vcf), value = 'nucleotide')
full_dict = dict_generator(all_cases,len(all_cases),value = 'nucleotide')


def index_annot(chr_haplotype,max_ind,chr_dictionary):

    # chr_haplotype: h*chr* matrix generated previosly
    # max_ind: Max of individuals in the vcf sample
    # chr_dictionary: chr*_dict variable to use

    temporal_indexes = [[]for x in xrange(0,max_ind)]

    for dlist_index,sublists in enumerate(chr_haplotype.tolist(), start= 0):
        for element_index,i in enumerate(sublists, start=0):
            if i == 1:
                temporal_indexes[dlist_index].append(element_index)# A way to
                # store indexes, counts every element inside a list and annotate
                # the count when it finds a 1.
    # Above iteration returns a list of lists with indexes were a match was
    # found
            elif i == 0:
                temporal_indexes[dlist_index].append('-')

    translated_indexes = [[]for x in xrange(0,max_ind)]

    for list_index, lists1 in enumerate(temporal_indexes, start= 0):
        if not lists1:
            pass
        else:
            for i in lists1:
                current_element = chr_dictionary[i]
                translated_indexes[list_index].append(current_element)

    return(translated_indexes)

    # Above iteration returns rsID/nucleotide translated version of the list of
    # lists generated in the
    # previous iteration.

# Made to manage the presence of '-' in dictionaries
# Returns a tuple containing np.arrays with unique string haplotype
# configurations and counts for its ocurrences
def almost_haplotypes(haplotype,max_ind,dictionary):
    annotate = index_annot(haplotype,max_ind,dictionary)

    temp_list = [[] for x in xrange(0,max_ind)]

    for num_list,ind_list in enumerate(annotate, start = 0):
        temp_list[num_list].append(''.join(ind_list))

    temp_list2 = np.unique(np.array(temp_list), return_counts = True)

    return(temp_list2)

# Final report at the stdout for partial haplotypes
def partial_haplotypes():

    # For haplotype 1
    trans_chr3st = np.unique((np.array(index_annot(h1chr3st,503,chr3st_dict))),
    return_counts = True)
    trans_chr3nd = np.unique((np.array(index_annot(h1chr3nd,503,chr3nd_dict))),
    return_counts = True)
    trans_chr5 = np.unique((np.array(index_annot(h1chr5,503,chr5_dict))),
    return_counts = True)
    trans_chr7 = np.unique((np.array(index_annot(h1chr7,503,chr7_dict))),
    return_counts = True)
    trans_chr20 = np.unique((np.array(index_annot(h1chr20,503,chr20_dict))),
    return_counts = True)

    # For haplotype 2

    strans_chr3st = np.unique((np.array(index_annot(h2chr3st,503,chr3st_dict))),
    return_counts = True)
    strans_chr3nd = np.unique((np.array(index_annot(h2chr3nd,503,chr3nd_dict))),
    return_counts = True)
    strans_chr5 = np.unique((np.array(index_annot(h2chr5,503,chr5_dict))),
    return_counts = True)
    strans_chr7 = np.unique((np.array(index_annot(h2chr7,503,chr7_dict))),
    return_counts = True)
    strans_chr20 = np.unique((np.array(index_annot(h2chr20,503,chr20_dict))),
    return_counts = True)

    print """These are the pieces of haplotype conserved for the TERC region
    >>> chromosome 3 <<<
    >>> string 1 <<<"""
    print trans_chr3st

    print """*These are the pieces of haplotype conserved for the TERC region*
    >>> chromosome 3 <<<
    >>> string 2 <<<"""
    print strans_chr3st

    print """*These are the pieces of haplotype conserved for the SENP7 region*
    >>> chromosome 3 <<<
    >>> string 1 <<<"""
    print trans_chr3nd

    print """*These are the pieces of haplotype conserved for the SENP7 region*
    >>> chromosome 3 <<<
    >>> string 2 <<<"""
    print strans_chr3nd

    print """*These are the pieces of haplotype conserved for the TERT region*
    >>> chromosome 5 <<<
    >>> string 1 <<<"""
    print trans_chr5

    print """*These are the pieces of haplotype conserved for the TERT region*
    >>> chromosome 5 <<<
    >>> string 2 <<<"""
    print strans_chr5

    print """*These are the pieces of haplotype conserved for the GPR37/POT1 region*
    >>> chromosome 7 <<<
    >>> string 1 <<<"""
    print trans_chr7

    print """*These are the pieces of haplotype conserved for the GPR37/POT1 region*
    >>> chromosome 7 <<<
    >>> string 2 <<<"""
    print strans_chr7

    print """*These are the pieces of haplotype conserved for the RTEL1 region*
    >>> chromosome 20 <<<
    >>> string 1 <<<"""
    print trans_chr20

    print """*These are the pieces of haplotype conserved for the RTEL1 region*
    >>> chromosome 20 <<<
    >>> string 2 <<<"""
    print strans_chr20


# Final report at the stdout for full haplotypes
def full_haplotypes():
    # For haplotype 1
    trans_full = almost_haplotypes(haplotype_1,503,full_dict)

    # For haplotype 2
    strans_full = almost_haplotypes(haplotype_2,503,full_dict)


    print """*This is the first haplotype string incluiding all the interest regions
    in all the case individuals*"""
    print trans_full

    print """This is the second haplotype string incluiding all the interest regions
    in all the case individuals*"""
    print strans_full


full_haplotypes()
