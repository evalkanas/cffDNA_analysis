#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
# Modified by Elise Valkanas <evalkanas@g.harvard.edu>

"""

"""


import argparse
import sys
import itertools
import pysam
from pysam import VariantFile
import re
import pdb

##def get_denovo_candidates(record, max_parents=10):
def get_all_denovo_candidates(record, proband, pat, mat):
    """
    Obtain list of samples which are putatively called de novo
    """
    #get list of genotypes for samples in the vcf
    #gts = [s['GT'] for s in record.samples.values()]
    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]    

    null_gt = [(None, None), (None,)]
    denovo = []
    #ignore the site if any individuals in the family have null genotypes at that position
    if len([i for i in null_gt if i in gts]) == 0:
    #fam vcf files always list individuals in the same order: proband, father, mother
    #find the possible mendelian inherited gentypes for the proband
#        mend_gts =  [ (gts[2][0],gts[1][0]), (gts[2][1],gts[1][0]), (gts[2][0],gts[1][1]), (gts[2][1],gts[1][1]) ] 
#        pot_gts = list(set([tuple(sorted(s)) for s in mend_gts])) # sort the order of each of these possible genotypes

    #check to see if the gt of the proband is one of these 
#        if gts[0] not in pot_gts:
        if gts[0] == (0,1) and gts[1] == (0,0) and gts[2] == (0,0):

            denovo.append(proband)

    return denovo


def get_male_Y_dn(record, proband, pat, mat):
    """
    Obtain list of samples which are putatively called de novo - uses known position of family members in vcf, but can easily use the IDs
    Do not call cases where the son's genotype is the same as the fathers on chrY - this is not de novo
    """
    #get list of genotypes for samples in the vcf
    #gts = [s['GT'] for s in record.samples.values()]
    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]

    null_gt = [(None, None), (None,)]
    denovo = []
    #ignore the site if the father or male child have null genotypes at that position
    if len([i for i in null_gt if i in gts]) == 0:
        mend_gts =  [ (gts[2][0],gts[1][0]), (gts[2][1],gts[1][0]), (gts[2][0],gts[1][1]), (gts[2][1],gts[1][1]) ]
        pot_gts = list(set([tuple(sorted(s)) for s in mend_gts]))
        if gts[0] not in pot_gts and gts[0] != gts[1]:
            denovo.append(proband)

    elif gts[2] in null_gt and gts[1] not in null_gt and gts[0] not in null_gt: #if only mother has null gt, we can still call this dn if son doesnt have same gt as dad
    #check to see if the gt of the proband is one of these
        if gts[0] != gts[1]:
            denovo.append(proband)
    return denovo

def get_male_X_dn(record, proband, pat, mat):
    """
    Obtain list of samples which are putatively called de novo - uses known position of family members in vcf, but can easily use the IDs
    Do not call cases where the son's genotype is two of the maternal alleles on chrX - this is not de novo
    """
    #get list of genotypes for samples in the vcf
    #gts = [s['GT'] for s in record.samples.values()]
    gts=[record.samples[proband]['GT'], record.samples[mat]['GT']]

    null_gt = [(None, None), (None,)]
    denovo = []
    #ignore the site if any individuals in the family have null genotypes at that position
    if len([i for i in null_gt if i in gts]) == 0:
    #fam vcf files always list individuals in the same order: proband, father, mother
    #find the possible mendelian inherited gentypes for the proband
#        mend_gts =  [ (gts[2][0],gts[1][0]), (gts[2][1],gts[1][0]), (gts[2][0],gts[1][1]), (gts[2][1],gts[1][1]) ]
#        pot_gts = list(set([tuple(sorted(s)) for s in mend_gts])) # sort the order of each of these possible genotypes

    #check to see if the gt of the proband is one of these
#        if gts[0] not in pot_gts and gts[0] != (gts[2][0],gts[2][0]) and  gts[0] != (gts[2][1],gts[2][1]):
        if gts[0] != (0,0) and gts[1] == (0,0):
            denovo.append(proband)
    return denovo

def get_hom_alt(record, proband, pat, mat, sex):
    """
    return if the variant is private to the family
    """
    #get list of genotypes for samples in the vcf
#    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]
    null_gt = [(None, None), (None,)]
    hom_alt = []
    #find positions where the proband is homalt and neither parent is hom alt
#    if gts[0] == (1,1) and (gts[1] != (1,1) and gts[2] != (1,1)): 
    if record.contig.lstrip("chr") == 'X' and sex == 1: # men on chr X
        gts = [record.samples[proband]['GT'], record.samples[mat]['GT']]
        if gts[0] == (1,1) and gts[1] == (0,1) and len([i for i in null_gt if i in gts]) == 0: # if male proband is hom and mom is het this is hom rec variant
            hom_alt.append(proband)
    elif record.contig.lstrip("chr") == 'X' and sex == 2: # females on chr X
        gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]
        if gts[0] == (1,1) and gts[1] != (0,0) and gts[2] == (0,1) and len([i for i in null_gt if i in gts]) == 0: # if female proband is hom alt and father is either het or hom alt and mother is het then this is a hom rec variant
            hom_alt.append(proband)
    else: # autosomes (and chrY) 
        gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]
        if gts[0] == (1,1) and (gts[1] == (0,1) and gts[2] == (0,1)) and len([i for i in null_gt if i in gts]) == 0:  # if proband is hom alt and both parents are het
            hom_alt.append(proband)
    return hom_alt

def get_proband_hets(record, proband, pat, mat):
    """
    return number of probands het at this site
    """
    null_gt = [(None, None), (None,)]
    het_gt = [(0, 1)]
    hets=[]
    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]
    # if proband is het and the parents are not no call ./. this site is a het
    if gts[0] in het_gt and len([i for i in null_gt if i in gts]) == 0: 
        hets.append(proband)
    return hets

def get_all_private(record, proband, pat, mat):
    """
    return if the variant is private to the family
    """
    #get list of genotypes for samples in the vcf
#    gts = [s['GT'] for s in record.samples.values()]
    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]

    null_gt = [(None, None), (None,)]
    private = []
    #ignore the site if any individuals in the family have null genotypes at that position
    if len([i for i in null_gt if i in gts]) == 0:
    #fam vcf files always list individuals in the same order: proband, father, mother
    #find positions where the proband is variant and has inherited this from only one parent
        if gts[0] == (0,1) and ((gts[1] == (0,1)) ^ (gts[2] == (0,1))):
            private.append(proband)
    return private

def get_all_private_x(record, proband, pat, mat, sex):
    """
    return private variants on the sex chromosomes
    """
    #get list of genotypes for samples in the vcf
#    gts = [s['GT'] for s in record.samples.values()]
    gts=[record.samples[proband]['GT'], record.samples[pat]['GT'], record.samples[mat]['GT']]

    null_gt = [(None, None), (None,)]
    private = []
    #ignore the site if any individuals in the family have null genotypes at that position
    if len([i for i in null_gt if i in gts]) == 0:
    #fam vcf files always list individuals in the same order: proband, father, mother
    #find positions where the proband is variant and has inherited this from only one parent
        if sex == 2:
            if gts[0] != (0,0) and ( ((gts[1] == (0,1)) ^ (gts[2] == (0,1))) or (gts[1] == (1,1) and gts[2]==(0,0)) ):
                private.append(proband)
        elif sex == 1:
            if gts[0] != (0,0) and ( (gts[1] == (0,1)) ^ (gts[2] == (0,1)) ): # only include spots where male proband is het (PAR) since all 1/1 (hemizygous) variants will be returned by the homozygous filter
                private.append(proband)
        else:
            print ("**** ERROR: Unrecognized sex: {} for proband: {} ****".format(sex, proband))
            exit(1)
    return private


def get_parent_vars(record, parent_list):
    """
    return number of parents that have this variant, variable will be zero if only parents have this variant
    """
    
    # not interested in individuals with ref or not called GT
    null_gt = [(0, 0), (None, None), (0, ), (None, )]

    #get list of probands with variant genotypes in the vcf
    parent_vars = [s for s in record.samples.keys() if record.samples[s]['GT'] not in null_gt and s in parent_list]
    types=[record.samples[s]['GT'] for s in parent_vars]
    AN=sum([int(list(s)[0])+int(list(s)[1]) for s in types]) # get allele number for parents with this variant by adding genotypes alleles (will only work for biallelic vars)

    #return the number of parentss variant at this position
    return AN

def get_proband_vars(record, proband_list):
    """
    return number of probands that have this variant, variable will be zero if only parents have this variant
    """

    # not interested in individuals with ref or not called GT
    null_gt = [(0, 0), (None, None), (0, ), (None, )]

    #get list of probands with variant genotypes in the vcf and add genotypes to list
    proband_vars = [record.samples[s]['GT'] for s in record.samples.keys() if record.samples[s]['GT'] not in null_gt and s in proband_list]

    #return the number of probands variant at this position
    return proband_vars

def get_ped_structure(ped_file): #make a dictionary for each individual of their parents (as long as both parents are not missing)
    ped_dict = {}
    founders = []
    affecteds = []
    with open(ped_file) as f:
        for line in f:
            fam, indiv, pat, mat, sex, aff = line.strip().split("\t")
            if aff == "2":
                affecteds.append(indiv)
#            if pat != "0" and  mat != "0": #consider only individuals with both parents in ped **Does not take affected status into account**
            ped_dict[indiv] = [pat, mat, sex] # create child dictionary
#                founders.extend([pat, mat])	# for all children with both parents, add those parents to a list
    return(ped_dict, founders, affecteds)

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='input vcf file - must be annotated')
    parser.add_argument('fout', help='output directory - to write dominant and homozygous variant files')
    parser.add_argument('proband', help='proband ID as it appears in the VCF and ped files') 
    parser.add_argument('ped', help='ped file')
    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin, 'r')
        name=sys.stdin.strip().split("/")[-1]
    else:
        vcf = pysam.VariantFile(args.vcf, 'r')
        name=args.vcf.strip().split("/")[-1]

    if args.fout in '- stdout'.split(): #just need 2 output files - recessive and het variants
#        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
        hom_file = pysam.VariantFile(sys.stdout+"/hom_variants_"+name, 'w', header=vcf.header)
        het_file = pysam.VariantFile(sys.stdout+"/het_variants_"+name, 'w', header=vcf.header)
        comp_het_file = pysam.VariantFile(sys.stdout+"/comp_het_variants_"+name, 'w', header=vcf.header)

    else:
#        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)
        hom_file = pysam.VariantFile(args.fout+"/hom_variants_"+name, 'w', header=vcf.header)
        het_file = pysam.VariantFile(args.fout+"/het_variants_"+name, 'w', header=vcf.header)
        comp_het_file = pysam.VariantFile(args.fout+"/comp_het_variants_"+name, 'w', header=vcf.header)

    ped_dict, parent_list, affected_list = get_ped_structure(args.ped) #get a dict of prodband: dad, mom, sex and a list of all parents of trio probands, and all affected individuals
    proband = args.proband
    pat = ped_dict[proband][0]
    mat = ped_dict[proband][1]
    sex = int(ped_dict[proband][2])

## TODO: add code to separate the trio probands from single probands and make singletons go through genotype filtering instead of inheritance filtering - have general filtering come first in script. 

    cohort_gt_cutoff = 4
    comp_het_list = [] 

 ### Set Variants ###
    cons_keep=['frameshift_deletion', 'frameshift_insertion', 'nonsynonymous_SNV', 'stopgain', 'stoploss', 'splicing'] #, 'UTR3', 'UTR5', 'upstream']
    cons_keep.extend(['nonframeshift_deletion', 'nonframeshift_insertion'])
    clin_remove=['Affects', 'Likely_benign', 'other', 'Benign/Likely_benign', 'drug_response', 'protective', 'association', 'risk_factor', '_protective', '_risk_factor', 'Benign', '_other', '_drug_response','_association']
    clin_keep=['Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic']

    for record in vcf:

        # Remove multiallelics
        if len(record.alts) > 1:
            continue

        # Do not include chr Y or MT 
        chrom = record.contig.lstrip("chr")

        if chrom  == 'MT' or chrom == 'Y':
            continue


        # remove spanning deletions
        if '*' in record.alts:
            continue

        if record.samples[proband]['GT'] == (1,1): #we are not interested in hom sites in mom since they are unlikley to cause disease
            continue


        # remove sites where the disease gene list annotation is not the same as the gene overlapping the consequence reported by annovar
        # TODO remove need to have disease_list_gene annotation and allow us to supply list of gene names for filtering
        disease_list=list(record.info['disease_list_gene'])
        gene_effect = list(record.info['AAChange.ensGene']) + list(record.info['GeneDetail.ensGene']) + list(record.info['Gene.ensGene'])
        ensGene = [s.strip().split(":")[0] for s in gene_effect]
        if len(set(disease_list).intersection(set(ensGene))) < 1:
            continue

        # get number of probands that are het and hom at this site
        candidates = get_proband_vars(record, affected_list)
        het_candidates = len([geno for geno in candidates if geno == (0,1)]) # number of probands that are het at this site
        hom_candidates = len([geno for geno in candidates if geno == (1,1)]) # numner of probands that are hom at this site

 ### Filter Variants ###

        ## Set variant and Gene list requirements
        # Is the variant of a consequence of interest according to gencode annotations
        cons=list(record.info['Func.ensGene']) + list(record.info['ExonicFunc.ensGene'])
        num_coding=len([e for e in cons if e in cons_keep])
        
        ## Allele frequency filtering
        dom_af_thresh = 0.01
        rec_af_thresh = 0.05
        raw_af_list = [record.info['ExAC_ALL'], record.info['AF_popmax'][0], record.info['1000g2015aug_all']] #Exac, gnomad, 1Kg
        raw_af_list.append(record.info['AF'][0]) #do we want to set same AF filters for cohort as dbs?

        # ANNOVAR doesn't give annotations for each multiallelic variant allele
        af_list = [ af for af in raw_af_list if (af is not None ) and ( af != "." ) ] #remove null values
        if not af_list: # if there are af values in the above 4 databases
            db_max_af = 0
        else:
            db_max_af = max(float(value) for value in af_list)
        max_af = db_max_af
        if max_af > rec_af_thresh: #af_threshold: # remove vars that are too frequent
            continue

        # remove if not coding effect of interest 
        if num_coding < 1 or record.info['disease_list_source'][0] == "NA": #remove if not in annotated gene list`:
            continue
        #if len(record.alts[0]) > 50 or len(record.ref) > 50: # remove indels larger than 50 bp
        #    continue

        # remove missense variants that are benign in clinvar
        if list(record.info['ExonicFunc.ensGene'])[0] == 'nonsynonymous_SNV' and set(record.info['CLNSIG']).issubset(clin_remove): #get rid of missense vars in clinvar as benign/risk factor, etc
            continue
        proband_format = dict(zip(record.format.keys(), record.samples[proband].values()))
        proband_alt_index = proband_format['GT'][1]
        
        if sum(list(proband_format['AD'])) == 0: #skip sites with no reads 
            continue

        AB=float(list(proband_format['AD'])[1])/float(sum(list(proband_format['AD']))) # calculate AB using the depth of the second allele in GT field (alt allele) 

        #### INHERITANCE FILTERING ####
        proband_gt = record.samples[proband]['GT']
        dom_list = []
        rec_list= []
#        print("inheritance: {}  AB: {} max_af: {}".format(record.info['disease_list_inheritance'], AB, max_af))
        if 'Dominant' in record.info['disease_list_inheritance'] and 0 < AB < 0.25 and max_af <= dom_af_thresh: 
            het_file.write(record)
            dom_list.append(record) 

        if 'Recessive' in record.info['disease_list_inheritance']:
            if 0 < AB < 0.75: 
                comp_het_list.append(record)
                if AB > 0.5:
                    rec_list.append(record)
                    hom_file.write(record)
        
    #remove half hets and print comp het list to file 
    #make list of all genes with hets that occur more than once in comp het list
    genes = []
    print([s.pos for s in comp_het_list])
    for s in comp_het_list:
        genes+=list(s.info['disease_list_gene'])
    # remove genes with only one occurance 
    keep_genes = list(set([i for i in genes if genes.count(i) > 1]))
    #loop through comp het list and keep only records where the gene annotation is in our non-unique genes list
    for record in comp_het_list:
        counter=0
        for gene_name in list(record.info['disease_list_gene']): #for each gene annotation, check to see if it is in our gene list and if so add to counter
            if gene_name in keep_genes:
                counter+=1
        if counter > 0: #if counter is greater than 0, at least one gene annotation is in our list to keep so write the record to the comp het list 
            comp_het_file.write(record)


    het_file.close()
    hom_file.close()
    comp_het_file.close()

if __name__ == '__main__':
    main()
