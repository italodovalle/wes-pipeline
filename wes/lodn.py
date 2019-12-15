#! /usr/bin/env python

#######################################################################
#		     Italo Faria do Valle                             #
#								      #
#		         lodn.py				      #
#		   italodovalle@gmail.com			      #
#######################################################################

import vcf
import pysam
import math
import pandas as pd
import re

def get_covs(samfile,pos,ch,alt,ref,min_mapq = 20, min_phred = 20):
    alt_cov = 0
    total_cov = 0
    ref_cov = 0
    ### A pileup column. A pileup column contains all the reads
    ### that map to a certain target base
    for pileupcolumn in samfile.pileup(str(ch),pos-2,pos):
        ####pileupcolumn.pos is 0-based
        if pileupcolumn.pos == pos - 1:      
            total_cov = pileupcolumn.n
            ###each read aligned to a pileup colum
            for pileupread in pileupcolumn.pileups:    
                b = pileupread.alignment.seq[pileupread.query_position]
                qvalue = pileupread.alignment.qual[pileupread.query_position]
                qvalue = ord(qvalue) - 33
                if pileupread.alignment.mapq >= min_mapq:
                    if qvalue >= min_phred:
                        if b == alt:
                            alt_cov += 1
                        elif b == ref:
                            ref_cov += 1
    return(alt_cov, ref_cov, total_cov)


def log_l_m_f(samfile,pos,ch,alt,ref,f,min_mapq = 20, min_phred = 20):
    prob = []
    pos = int(pos)
    for pileupcolumn in samfile.pileup(str(ch),pos - 2,pos):
        if pileupcolumn.pos == pos - 1 :
            for pileupread in pileupcolumn.pileups:
                b = pileupread.alignment.seq[pileupread.query_position]
                qvalue = pileupread.alignment.qual[pileupread.query_position]
                qvalue = ord(qvalue) - 33
                if pileupread.alignment.mapq >= min_mapq:
                    if qvalue >= min_phred:
                        e = 10**(float(-1 * qvalue)/10)
                        if b == ref:
                            p = (f*(e/3)) + ((1-f)*(1-e))
                        elif b == alt:
                            p = (f*(1-e) )+((1-f)*(e/3))
                        else:
                            p = e/3
                        prob.append(math.log10(p))
    if len(prob) == 0:
        return (float('NaN'))
    else:
        return (sum(prob))


def lodn_calculator_vcf(tumor_vcf,normal_bam,outfile,min_n_cov = 8,
                        min_t_cov = 14,lodn_known = 5.5, lodn_novel=2.2):  
    bamfile = pysam.Samfile(normal_bam, "rb")
    vcf_reader = vcf.Reader(open(tumor_vcf, 'r'))
    vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader)
    
    for variant in vcf_reader:
        if variant.is_snp and len(variant.ALT) == 1:
            tumor_cov = variant.INFO['DP']
            #if not variant.ID:
            #    dbsnp = 'NOVEL'
            l_m0 = 0
            l_m = 0
            alt_cov, ref_cov, total_cov = get_covs(bamfile,variant.POS,
                                                   variant.CHROM,
                                                   variant.ALT[0],
                                                   variant.REF)
            if total_cov > min_n_cov and tumor_cov > min_t_cov:
                l_m0 = log_l_m_f(bamfile,variant.POS,variant.CHROM,
                                 variant.ALT[0],variant.REF,0)
                l_m = log_l_m_f(bamfile,variant.POS,variant.CHROM,
                                variant.ALT[0],variant.REF,0.5)
                try:
                    lod_n = (l_m0)-(l_m)
                except:
                    lod_n = 0
                if not variant.ID:
                    if lod_n >= lodn_known:
                        vcf_writer.write_record(variant)
                else:
                    if lod_n >= lodn_novel:
                        vcf_writer.write_record(variant)

def lodn_calculator_table(table, sep, normal_bam ,outfile=None, min_n_cov = 8, 
                          min_t_cov = 14,lodn_known = 5.5, lodn_novel=2.2):
    
    patient = pd.read_table(table,sep=sep)
    samfile = pysam.Samfile(normal_bam, "rb")
    ix = []
    ## creating columns
    patient['t_cov.ref'] = 0; patient['t_cov.alt']=0
    patient['n_cov.ref']=0; patient['n_cov.alt']=0
    patient['filter'] = 'LODn'
    
    z = []
    for label in patient.index:
        alt = patient['alt.base'].loc[label]
        ref = patient['ref.base'].loc[label]
        ch = patient['chr'].loc[label]
        pos = int(patient['pos.start'].loc[label])
        tumor_cov = patient['cov.ref'].loc[label] + patient['cov.alt'].loc[label]
        patient['t_cov.ref'].loc[label] = patient['cov.ref'].loc[label]
        patient['t_cov.alt'].loc[label] = patient['cov.alt'].loc[label]
        patient['lodn'] = 0
        l_m0 = 0
        l_m = 0
        ref_cov = 0
        ## Check if it is a single nt variant
        if len(alt) == 1 and len(ref) == 1:
            if re.match('[A,C,T,G]',alt) and re.match('[A,C,T,G]',ref):
                alt_cov, ref_cov, total_cov = get_covs(samfile,pos,ch,alt,ref)
                if total_cov > min_n_cov and tumor_cov > min_t_cov:
                    patient['n_cov.ref'].loc[label] = ref_cov
                    patient['n_cov.alt'].loc[label] = alt_cov
                    l_m0 = log_l_m_f(samfile,pos,ch, alt, ref, 0)
                    l_m = log_l_m_f(samfile, pos, ch, alt, ref, 0.5)
                try:
                    lod_n = (l_m0)-(l_m)
                except:
                    lod_n = 0            
                patient['lodn'].loc[label] = lod_n
                if patient['known.flag'].loc[label] == 1:      
                    if lod_n >= lodn_known:
                        ix.append(label)
                elif patient['known.flag'].loc[label] == 0:
                    if lod_n >= lodn_novel:
                        ix.append(label)
       
    patient = patient.loc[ix]
    if outfile:
        patient.to_csv(outfile, sep=sep, header=True, index=True)
    
    return (patient)
