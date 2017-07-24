import pandas as pd
import numpy as np
import os
import sys

# This is a final version of the haploseparator, made to run in cluster
# managing ~4000 individuals VCF file.

# NO 1K GENOMES REPORT FOR THIS VERSION

#Version made to work with the hrc derived haplotypes, directory1 variable changed!

# GENERAL FUCNTIONS COLLECTION

global max_individuals_1k, max_individuals_hrc_cases, max_individuals_hrc_controls
global directory1
max_individuals_1k = 503
max_individuals_hrc_cases = 4328
max_individuals_hrc_controls = 7046
directory1 = '/panfs/pan5/drobles/igarcia/'

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

# Next function is made to manage the presence of '-' in dictionaries
# Returns a tuple containing np.arrays with unique string haplotype
# configurations and counts for its ocurrences
def almost_haplotypes(haplotype,max_ind,dictionary):
    annotate = index_annot(haplotype,max_ind,dictionary)

    temp_list = [[] for x in xrange(0,max_ind)]

    for num_list,ind_list in enumerate(annotate, start = 0):
        temp_list[num_list].append(''.join(ind_list))

    raw_report = np.unique(np.array(temp_list), return_counts = True)

    most_common_haplotypes = raw_report[0].tolist()
    most_common_counts = raw_report[1].tolist()
    common_haplotypes = sorted(zip(most_common_haplotypes,most_common_counts),
    key=lambda seq_string: seq_string[1],reverse=True)[:3]

    return(common_haplotypes)

def each_chromosome():
    hrc_vcfs_list_ca = []
    hrc_vcfs_list_co = []
    regions = ['RTEL1_case_derived', 'RTEL1_control_derived']
    for hrc_vcfs_ca in os.listdir(directory1):
        if hrc_vcfs_ca[-3:] == 'vcf' and hrc_vcfs_ca[:6] == 'hrc_ca':
            hrc_vcfs_list_ca.append(hrc_vcfs_ca)
    for hrc_vcfs_co in os.listdir(directory1):
        if hrc_vcfs_co[-3:] == 'vcf' and hrc_vcfs_co[:6] == 'hrc_co':
            hrc_vcfs_list_co.append(hrc_vcfs_co)
    for vcf_name_ca,vcf_name_co,re_names in zip(hrc_vcfs_list_ca,hrc_vcfs_list_co,regions):

        pandas_df_cas = pd.read_table(directory1+vcf_name_ca, sep = '\t')
        pandas_df_con = pd.read_table(directory1+vcf_name_co, sep = '\t')

        current_dict = dict_generator(pandas_df_cas,len(pandas_df_cas),
        value = 'nucleotide')


        haplo1_cas = haploseparator(pandas_df_cas,len(pandas_df_cas),
        haplotype = 'h1')
        haplo2_cas = haploseparator(pandas_df_cas,len(pandas_df_cas),
        haplotype = 'h2')

        haplo1_con = haploseparator(pandas_df_con,len(pandas_df_con),
        haplotype = 'h1')
        haplo2_con = haploseparator(pandas_df_con,len(pandas_df_con),
        haplotype = 'h2')

        translated_haplo1_cas = almost_haplotypes(haplo1_cas,
        max_individuals_hrc_cases,current_dict)
        translated_haplo2_cas = almost_haplotypes(haplo2_cas,
        max_individuals_hrc_cases,current_dict)

        translated_haplo1_con = almost_haplotypes(haplo1_con,
        max_individuals_hrc_controls,current_dict)
        translated_haplo2_con = almost_haplotypes(haplo2_con,
        max_individuals_hrc_controls,current_dict)

        print """
        *These are the most common haplotype configurations in the
        first string of haplotype for the """,re_names,""" region*
        >>>CASE INDIVIDUALS<<<
        """, translated_haplo1_cas
        print """
        *These are the most common haplotype configurations in the
        second string of haplotype for the """,re_names,""" region*
        >>>CASE INDIVIDUALS<<<
        """, translated_haplo2_cas
        print """
        *These are the most common haplotype configurations in the
        first string of haplotype for the """,re_names,""" region*
        >>>CONTROL INDIVIDUALS<<<
        """, translated_haplo1_con
        print """
        *These are the most common haplotype configurations in the
        second string of haplotype for the """,re_names,""" region*
        >>>CONTROL INDIVIDUALS<<<
        """, translated_haplo2_con

    to_full= raw_input("Include full haplotype report?(y/n): ")
    if to_full == 'y':
        all_chromosomes_dis()
    if to_full == 'n':
        print ("*END*")

def all_chromosomes_dis():
    #FULL VCF FOR ALL REGIONS IN CASES
    ##NOT USE FOR TERC TEST!
    all_cases = pd.read_table(directory1+'full_diseased_w_o_header.vcf',
    sep='\t')
    all_controls = pd.read_table(directory1+'full_controls_w_o_header.vcf',
    sep='\t')
    haplotype_1 = haploseparator(all_cases,len(all_cases),
    haplotype = 'h1')
    haplotype_2 = haploseparator(all_cases,len(all_cases),
    haplotype = 'h2')
    haplotype_1c = haploseparator(all_controls,len(all_controls),
    haplotype = 'h1')
    haplotype_2c = haploseparator(all_controls,len(all_controls),
    haplotype = 'h2')
    full_dict = dict_generator(all_cases,len(all_cases),
    value = 'nucleotide')

    #For haplotype 1
    trans_full = almost_haplotypes(haplotype_1,
    max_individuals_hrc_cases,full_dict)

    print """
*These are the most common haplotype configurations in the
first string of haplotype for all the interest regions*
>>>CASE INDIVIDUALS<<<
""", trans_full

    #For haplotype 2
    strans_full = almost_haplotypes(haplotype_2,
    max_individuals_hrc_cases,full_dict)

    print """
*These are the most common haplotype configurations in the
first string of haplotype for all the interest regions*
>>>CASE INDIVIDUALS<<<
""", strans_full

# The same just for controls

    #For haplotype 1 (controls)
    trans_fullc = almost_haplotypes(haplotype_1c,
    max_individuals_hrc_controls,full_dict)

    print """
*These are the most common haplotype configurations in the
first string of haplotype for all the interest regions*
>>>CONTROL INDIVIDUALS<<<
""", trans_fullc

    #For haplotype 2 (controls)
    strans_fullc = almost_haplotypes(haplotype_2c,
    max_individuals_hrc_controls,full_dict)

    print """
*These are the most common haplotype configurations in the
first string of haplotype for all the interest regions*
>>>CONTROL INDIVIDUALS<<<
""", strans_fullc

    print ("*END*")

def helper():
    partialvsfull = raw_input('''
Show full haplotype report or partial haplotype report (per region)?(f/p): ''' )
    print "*Working with the 4328 cases from the HRC imputed melanoma GWAS*"
    print "*Working with the 7046 controls from the HRC imputed melanoma GWAS*"
    if partialvsfull == 'f':
        all_chromosomes_dis()
    if partialvsfull == 'p':
        each_chromosome()

helper()
