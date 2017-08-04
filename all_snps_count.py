import pandas as pd
import numpy as np
import os
import sys

test_vcf = pd.read_table("~/Escritorio/lead_snps_vcfs/haplotype_snps/cases_hrc_haplo_chr3st.vcf", sep = '\t')

# Preliminar variables

subs = test_vcf.iloc[:,9:] # Show only those columns and rows with haplotype info
individuals = subs.columns.values #Column names, one per individual
nindividuals = len(individuals) # Individuals in the present dataset
nsnps = (len(subs.index)-1) # SNPs in the region

# Next, lists of lists that appends a "string of haplotype" for each individual,
# sites of interest (rows).
haplotype1 = [[]for x in xrange(0,nindividuals)]
haplotype2 = [[]for x in xrange(0,nindividuals)]

# Separating the phased dataset into individual haplotype strings
# Appending the data in its respective ind list.
for index_ind,individual in enumerate(individuals, start = 0):
    for snp in range(0,nsnps):
        haplotype1[index_ind].append(subs[individual][snp][0])

for index_ind,individual in enumerate(individuals, start = 0):
    for snp in range(0,nsnps):
        haplotype2[index_ind].append(subs[individual][snp][2])

#Return the individual haplotype strings to a df format.

haplo1_df = pd.DataFrame(haplotype1)
haplo1_df = haplo1_df.transpose()
haplo1_df.columns = individuals

haplo2_df = pd.DataFrame(haplotype2)
haplo2_df = haplo2_df.transpose()
haplo2_df.columns = individuals

# One list for each site for all individuals in the dataset.
# One per haplotype.
h1_compilated_sites = [[] for x in xrange(0,nsnps)]
h2_compilated_sites = [[] for x in xrange(0,nsnps)]

# Indexing in second dataframe and annotating each site in a the previous lists.
for individual in individuals:
    for snp_index,snp in enumerate(range(0,nsnps), start = 0):
        h1_compilated_sites[snp_index].append(haplo1_df[individual][snp])

for individual in individuals:
    for snp_index,snp in enumerate(range(0,nsnps), start = 0):
        h2_compilated_sites[snp_index].append(haplo2_df[individual][snp])

# One list for each site, containing the number of individuals with:
# a) Reference allele b) Alternative allele and c) Deletions for the site.
h1_report = [[]for x in xrange(0,nsnps)]
h2_report = [[]for x in xrange(0,nsnps)]
homozygosity_report = [[]for x in xrange(0,nsnps)]

# Each list at the h*_compilated_sites has the allele value for all the
# individuals.
# At haplotype string 1
for index_site,site in enumerate(h1_compilated_sites, start = 0):
    REF_h1 = 0 # Restart the count at each site iteration
    ALT_h1 = 0
    DEL_h1 = 0
    for i in site: # each i represent an individual.
        if i == '0':
            REF_h1 += 1
        elif i == '1':
            ALT_h1 += 1
        else:
            DEL_h1 += 1
    h1_report[index_site].extend((REF_h1,ALT_h1,DEL_h1))

# At haplotype string 2
for index_site,site in enumerate(h2_compilated_sites, start = 0):
    REF_h2 = 0 # Restart the count at each site iteration
    ALT_h2 = 0
    DEL_h2 = 0
    for i in site: # each i represent an individual.
        if i == '0':
            REF_h2 += 1
        elif i == '1':
            ALT_h2 += 1
        else:
            DEL_h2 += 1
    h2_report[index_site].extend((REF_h2,ALT_h2,DEL_h2))

for index_site,site_1,site_2 in zip(range(0,nsnps),h1_compilated_sites,h2_compilated_sites):
    REF_homozygous = 0
    ALT_homozygous = 0
    heterozygous = 0
    for allele_a, allele_b in zip(site_1,site_2):
        if allele_a == '0' and allele_b == '0':
            REF_homozygous += 1
        elif allele_a == '0' and allele_b == '1':
            heterozygous += 1
        elif allele_a == '1' and allele_b == '1':
            ALT_homozygous += 1
        elif allele_a == '1' and allele_b == '0':
            heterozygous += 1
    homozygosity_report[index_site].extend((REF_homozygous,ALT_homozygous,heterozygous))

df_h1_report = pd.DataFrame(h1_report)
df_h2_report = pd.DataFrame(h2_report)
df_homozygosity_report = pd.DataFrame(homozygosity_report)
extra_fields = test_vcf[['ID','REF','ALT']]
pre_final_report = pd.concat([df_h1_report,df_h2_report,df_homozygosity_report], axis = 1)
final_report = pre_final_report.multiply(100).div(nindividuals)
out_report = pd.concat([extra_fields,final_report], axis = 1)
final_col_names = ['ID','REF','ALT','H1_REF_ALLELE','H1_ALT_ALLELE','H1_DEL','H2_REF_ALLELE',
'H2_ALT_ALLELE','H2_DEL','REF_HOMOZYGOUS','ALT_HOMOZYGOUS','HETEROZYGOUS']
out_report.columns = final_col_names

out_report

return(out_report)

hrc_vcfs_list_ca = []
hrc_vcfs_list_co = []
regions = ['RTEL1','OBCF1','TERC','TERT','GPR37/POT1']
for hrc_vcfs_ca in os.listdir(working_directory):
    if hrc_vcfs_ca[-3:] == 'vcf' and hrc_vcfs_ca[:3] == 'cas':
        hrc_vcfs_list_ca.append(hrc_vcfs_ca)
for hrc_vcfs_co in os.listdir(working_directory):
    if hrc_vcfs_co[-3:] == 'vcf' and hrc_vcfs_co[:3] == 'con':
        hrc_vcfs_list_co.append(hrc_vcfs_co)

for vcf_name_ca,vcf_name_co,re_names in zip(hrc_vcfs_list_ca,hrc_vcfs_list_co,regions):
    pandas_df_cas = pd.read_table(working_directory+vcf_name_ca, sep = '\t')
    pandas_df_con = pd.read_table(working_directory+vcf_name_co, sep = '\t')

    cases_report = haploreference(pandas_df_cas)
    controls_report = haploreference(pandas_df_con)
