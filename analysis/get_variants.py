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
#import svtk.utils as svu
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
            if pat != "0" and  mat != "0": #consider only individuals with both parents in ped **Does not take affected status into account**
                ped_dict[indiv] = [pat, mat, sex] # create child dictionary
                founders.extend([pat, mat])	# for all children with both parents, add those parents to a list
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

    if args.fout in '- stdout'.split():
#        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
        rec_homs = pysam.VariantFile(sys.stdout+"/man_rev_homs_"+name, 'w', header=vcf.header)
        manual_rev = pysam.VariantFile(sys.stdout+"/man_rev_dom_"+name, 'w', header=vcf.header)
        comp_het = pysam.VariantFile(sys.stdout+"/comp_het_"+name, 'w', header=vcf.header)
        private_vars = pysam.VariantFile(sys.stdout+"/private_vars_"+name, 'w', header=vcf.header)
    else:
#        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)
        rec_homs = pysam.VariantFile(args.fout+"/man_rev_homs_"+name, 'w', header=vcf.header)
        manual_rev = pysam.VariantFile(args.fout+"/man_rev_dom_"+name, 'w', header=vcf.header)
        comp_het = pysam.VariantFile(args.fout+"/comp_het_"+name, 'w', header=vcf.header)
        private_vars = pysam.VariantFile(args.fout+"/private_vars_"+name, 'w', header=vcf.header)

    ped_dict, parent_list, affected_list = get_ped_structure(args.ped) #get a dict of prodband: dad, mom, sex and a list of all parents of trio probands, and all affected individuals
    proband = args.proband
    print("proband {}".format(proband))
    #print(ped_dict)
    print("proband entry {}".format(ped_dict[proband]))
    pat = ped_dict[proband][0]
    mat = ped_dict[proband][1]
    sex = int(ped_dict[proband][2])

#    parent_list=[s[0] for s in list(ped_dict.values())]
#    parent_list.extend([s[1] for s in list(ped_dict.values())])

## TODO: add code to separate the trio probands from single probands and make singletons go through genotype filtering instead of inheritance filtering - have general filtering come first in script. 

    hom_num=0
    de_novo_num=0
    #cohort_gt_cutoff = 5
    cohort_gt_cutoff = 8

 ### Set Variants ###
    cons_keep=['frameshift_deletion', 'frameshift_insertion', 'nonsynonymous_SNV', 'stopgain', 'stoploss', 'splicing'] #, 'UTR3', 'UTR5', 'upstream']
    clin_remove=['Affects', 'Likely_benign', 'other', 'Benign/Likely_benign', 'drug_response', 'protective', 'association', 'risk_factor', '_protective', '_risk_factor', 'Benign', '_other', '_drug_response','_association']
    clin_keep=['Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic']

    for record in vcf:

        denovos = []
        homs = []
        hets =[] 
        private = []

        # Remove multiallelics
        if len(record.alts) > 1:
            continue

        # Do not include chr Y or MT (these shouldn't be called in the prenatals anyway)
#        if record.contig.strip('chr') == 'Y' or record.contig.strip('chr') == 'MT':
        if record.contig.strip('chr') == 'MT':
            continue

        # remove sites that overlap a deletion
        if '*' in record.alts:
            continue

        # remove no call sites in proband
        if None in record.samples[proband]['GT']:
            continue

        # remove locations where proband is hom ref
        if record.samples[proband]['GT'] == (0,0) or record.samples[proband]['GT'] == (None, ):
            continue

        # remove sites where the disease gene list annotation is not the same as the gene overlapping the consequence reported by annovar
        disease_list=list(record.info['disease_list_gene'])
        gene_effect = list(record.info['AAChange.ensGene']) + list(record.info['GeneDetail.ensGene']) + list(record.info['Gene.ensGene'])
        ensGene = [s.strip().split(":")[0] for s in gene_effect]
        if len(set(disease_list).intersection(set(ensGene))) < 1:
            continue

        # get number of probands that are het and hom at this site
        candidates = get_proband_vars(record, affected_list)
        het_candidates = len([geno for geno in candidates if geno == (0,1)]) # number of probands that are het at this site
        hom_candidates = len([geno for geno in candidates if geno == (1,1)]) # numner of probands that are hom at this site
#        if candidates > 5:
#            continue

        # remove sites where any family member has less than AD of 5
        ADs = [record.samples[proband]['AD'], record.samples[pat]['AD'], record.samples[mat]['AD']]
        pass_AD = [s for s in ADs if s != (None,) and sum(s) > 5]
        if len(pass_AD) < 3: # if all 3 individuals in this fam do not have passing AD, remove the site
            continue

        # remove sites in too many parents
        
        parent_vars = get_parent_vars(record, parent_list)
        if parent_vars > (len(parent_list)*2*0.05): # if more than 5% of parental alleles (across the cohort) are alt at this location, remove
#        if parent_vars > 332:    # if there are less than 5% alt allele in the parent pop - specific to ASD cohort of 3332 children
            continue


        #### INHERITANCE FILTERING ####
        # get de novo variants in this family - sex specific
        chrom = record.contig.lstrip("chr")

        if chrom  == 'MT' or chrom == 'Y':  # remove mito 
            continue

        if chrom == "X" and sex == 1:
            denovos = get_male_X_dn(record, proband, pat, mat)
        elif chrom  == "Y" and sex == 1:
            denovos = get_male_Y_dn(record, proband, pat, mat)
        else:
            denovos = get_all_denovo_candidates(record, proband, pat, mat)

        homs = get_hom_alt(record, proband, pat, mat, sex)
        hets = get_proband_hets(record, proband, pat, mat)
        
        # Get inherited variants in this family
        if chrom == 'X':
            private = get_all_private_x(record, proband, pat, mat, sex)
        else:
            # Skip records without any variants fitting genotype criteria
            private = get_all_private(record, proband, pat, mat)

        # Skip records without any variants of interest
        if len(denovos) == 0 and len(homs) == 0 and len(hets) == 0 and private == 0: #if the proaband isn't variant here, skip this record
            continue

 ### Filter Variants ###

        ## Set variant and Gene list requirements
        # Is the variant of a consequence of interest according to gencode annotations
        cons=list(record.info['Func.ensGene']) + list(record.info['ExonicFunc.ensGene'])
        num_coding=len([e for e in cons if e in cons_keep])
        
        ## Allele frequency filtering
#        af_threshold = 0.05
        dom_af_thresh = 0.01
        rec_af_thresh = 0.05
        raw_af_list = [record.info['ExAC_ALL'][0], record.info['gnomAD_genome_ALL'][0], record.info['gnomAD_exome_ALL'][0], record.info['1000g2015aug_all']]

        # ANNOVAR doesn't give annotations for each multiallelic variant allele
        af_list = [ af for af in raw_af_list if (af is not None ) and ( af != "." ) ] #remove null values
        if not af_list: # if there are af values in the above 4 databases
            db_max_af = 0
        else:
            db_max_af = max(float(value) for value in af_list)

#        max_af = max(db_max_af, record.info['AF_founder'])
        max_af = db_max_af

        if max_af > rec_af_thresh: #af_threshold: # remove vars that are too frequent
            continue

        # remove if not coding effect of interest 
        if num_coding < 1:   # or record.info['disease_list_source'][0] == "NA": remove if not in annotated gene list`
            continue

        # only keep PASS/ExcessHet SNVs
        if len(record.alts[0]) == 1 and len(record.ref) == 1:  #  SNV
            if record.filter.keys()[0] not in ['PASS', 'ExcessHet']: # only keep SNV vars that pass vqsr or were flagged as Excess hets - keep all indels regardless of filter
                continue
        else: # indel
            if len(record.alts[0]) > 50 or len(record.ref) > 50: # remove indels larger than 50 bp
                continue

        # remove missense variants that are benign in clinvar
        if list(record.info['ExonicFunc.ensGene'])[0] == 'nonsynonymous_SNV' and set(record.info['CLNSIG']).issubset(clin_remove): #get rid of missense vars in clinvar as benign/risk factor, etc
            continue

        proband_format = dict(zip(record.format.keys(), record.samples[proband].values()))
        proband_alt_index = proband_format['GT'][1]
#        if sum(list(proband_format['AD'])) < 5: # if there are less than 5 reads at this position, remove this variant      
#            continue

        AB=float(list(proband_format['AD'])[proband_alt_index])/float(sum(list(proband_format['AD']))) # calculate AB using the depth of the second allele in GT field / sum depths for both alleles so AB =1 for both 0/0 and 1/1 homs 
        if AB < 0.15: # remove sites with low AB - this will remove 0/0 de novo sites
            continue
        proband_gq=proband_format['GQ']
        if proband_gq is None: #or record.filter.keys()[0] != 'PASS':
            continue

#       convert CADD to float    
        CADD = record.info['CADD_phred'][0]
        if CADD == '.':
            CADD = 0
        CADD = float(CADD)


        if record.info['reg_mis_constr'] == 'mis_constrained_region':
            mis_constr = True
        else:
            mis_constr = False

#       missense filtering
        if list(record.info['ExonicFunc.ensGene'])[0] == 'nonsynonymous_SNV':
#       Tier 3 if missense variant P/LP in clinvar
            if len([s for s in list(record.info['CLNSIG']) if s in clin_keep]) > 0:
                mis_tier = 3
#           Tier 2 cadd > 30 or > 15 + missense constr
            elif CADD > 30 or (CADD > 15 and mis_constr):
                mis_tier = 2 
#           Tier 1 CADD > 15
            elif CADD > 15:
                mis_tier = 1
            else:
                mis_tier = 0 

        if len(private) > 0 and max_af <= dom_af_thresh and parent_vars <=  (len(parent_list)*2*0.01) and 'Dominant' in record.info['disease_list_inheritance']: #remove if no private variants at this site, if too common in ref databases or if in more than 1% pf parental alleles or if not in our annotated gene list
#        if len(private) > 0 and max_af <= dom_af_thresh and record.info['disease_list_source'][0] != "NA": #remove if no private variants at this site, if too common in ref databases or if in more than 1% pf parental alleles or if not in our annotated gene list
            if chrom == "X":
                if hom_candidates + het_candidates <= cohort_gt_cutoff and (list(record.info['ExonicFunc.ensGene'])[0] != 'nonsynonymous_SNV' or mis_tier >= 3):
                    private_vars.write(record)
            else:
                if het_candidates <= cohort_gt_cutoff and (list(record.info['ExonicFunc.ensGene'])[0] != 'nonsynonymous_SNV' or mis_tier >= 3):
                    private_vars.write(record) 

        if len(hets) > 0 and hom_candidates <= cohort_gt_cutoff and 'Recessive' in record.info['disease_list_inheritance']: # add site to compound het list if proband is het here and 5 or less probands are also het here
            comp_het.write(record)

        if len(denovos) > 0 and 'Dominant' in record.info['disease_list_inheritance']: # if variant is de novo and gene is associated with AD disease
            # GT specific filtering logic
            if int(proband_format['GT'][0]) + int(proband_format['GT'][1]) != 1 and hom_candidates <= cohort_gt_cutoff and (list(record.info['ExonicFunc.ensGene'])[0] != 'nonsynonymous_SNV' or mis_tier >= 1):  # if this variant is homozygous (ref or alt)  add to  file- all alleles here should be 1 and therefore all hets 0/1 and homs 1/1 or 0/0 so we can add these alleles to infer GT
                manual_rev.write(record)
                de_novo_num+=1
#                continue
            elif max_af <= dom_af_thresh and parent_vars <= (len(parent_list)*2*0.01) and het_candidates <= cohort_gt_cutoff and (list(record.info['ExonicFunc.ensGene'])[0] != 'nonsynonymous_SNV' or mis_tier >= 1): # if variant is het, apply stricter af filters before adding to file and remove if too common in ref databases, or if in more than 1% of parental alleles 
                manual_rev.write(record)
                de_novo_num+=1
#                continue
        if len(homs) > 0 and hom_candidates <= cohort_gt_cutoff and 'Recessive' in record.info['disease_list_inheritance'] and (list(record.info['ExonicFunc.ensGene'])[0] != 'nonsynonymous_SNV' or mis_tier >= 1): #for this record, if the individual has a hom variant here and less than or equal to 5 probands also have this variant
            rec_homs.write(record)
            hom_num+=1
#            continue
        else:
            continue

    manual_rev.close()
    rec_homs.close()
    comp_het.close()
    private_vars.close()
    print("Homs: {} De novo: {}".format(hom_num, de_novo_num))

if __name__ == '__main__':
    main()
