#! /usr/bin/env python

#######################################################################
#		     Italo Faria do Valle                             #
#								      #
#		        gatk_pipe.py				      #
#		   italodovalle@gmail.com			      #
#######################################################################

import os
import sys
import argparse
import pandas as pd
import yaml
from wes import ExomeFunctions as ef

if __name__ == '__main__':
    description = 'Whole Exome Analysis Pipeline'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', required=True, dest='fastq',nargs='*',
                        action='store',help='fastq files')
    parser.add_argument('-s', required=True, dest='name',action='store',
                        help='Sample name')
    parser.add_argument('-c', dest='cores', action='store',
                        help='number of cores to use, default = 2',
                        default= 2)
    parser.add_argument('-f', required=True, dest='config', action='store',
                        help='config file yaml')
    parser.add_argument('-d', dest='outdir', action='store',
                        help='output directory', default= '.')

    if len(sys.argv) <= 1: 
        parser.print_help() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    

    ##--------------------------------
    #Input file fastq already trimmed
    fastq = args.fastq
    name = args.name
    cores = int(args.cores)
    config = args.config
    outdir = os.path.abspath(args.outdir)
        

    with open(config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    sam = ef.Mapping(fastq,cfg['softwares']['bwa'],cfg['ref-files']['hg'],
                     outdir=outdir,name=name,
                     PL=cfg['sample-details']['platform'],
                     LB=cfg['sample-details']['library'],
                     cores=cores,err_log = True)
    
    bam = ef.Sort(sam,cfg['softwares']['picard'],tmp_dir=cfg['tmp_dir'],
                  outdir=outdir,name=name,err_log = True)
    
    marked = ef.MarkDuplicates(bam,cfg['softwares']['picard'],
                               tmp_dir=cfg['tmp_dir'],outdir=outdir,
                               name=name,err_log=True)
    
    ef.BuildBamIndex (marked,cfg['softwares']['picard'],
                      tmp_dir=cfg['tmp_dir'], outdir=outdir,
                      name=name, err_log=True)
    
    realigned = ef.IndelRealigner (marked, cfg['softwares']['gatk'],
                                   cfg['ref-files']['hg'],
                                   cfg['ref-files']['target_realigner'],
                                   cfg['ref-files']['indels_ref'],
                                   outdir=outdir,name=name,err_log=True,
                                   tmp_dir=cfg['tmp_dir'])
    
    recal = ef.BQSR(realigned, cfg['softwares']['gatk'], 
                    cfg['ref-files']['hg'],cfg['ref-files']['dbsnp'],
                    cfg['ref-files']['indels_ref'],name=name,cores=cores,
                    outdir=outdir,err_log=True,tmp_dir=cfg['tmp_dir'])
    
    genotyping_mode = cfg['hap-caller']['mode']
    stand_emit_conf=cfg['hap-caller']['emit_conf']
    stand_call_conf=cfg['hap-caller']['call_conf']
    vcf_raw = ef.HaplotypeCaller(recal,cfg['softwares']['gatk'],
                                 cfg['ref-files']['hg'],
                                 outdir=outdir,name=name,err_log=True,
                                 genotyping_mode = genotyping_mode,
                                 stand_emit_conf=stand_emit_conf,
                                 stand_call_conf=stand_call_conf)
    filter_exp_snps = cfg['hard-filter']['snps']
    filter_exp_indels = cfg['hard-filter']['indels']
    vcf_filt = ef.HardFilter (vcf_raw,cfg['softwares']['gatk'],
                              cfg['ref-files']['hg'],
                              filter_exp_snps, 
                              filter_exp_indels,
                              outdir=outdir,name=name,err_log=True)

    ann = ef.Annotation(vcf_filt,cfg['softwares']['convert2annovar'],
                        cfg['softwares']['annovar'],
                        cfg['ref-files']['humandb'],
                        cfg['ref-files']['build_ver'],
                        cfg['ref-files']['dbsnp_ver'],
                        cfg['ref-files']['kg_ver'],
                        cfg['ref-files']['mitochondrial_ver'],
                        outdir=outdir,name=name,err_log=True,maf=0.05)
    ef.MakeFinalFile (ann[0],ann[1],ann[2],outfile='%s/%s.tsv'%(outdir,name))
