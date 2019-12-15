#! /usr/bin/env python

#######################################################################
#		     Italo Faria do Valle                             #
#								      #
#		        mutect_pipe.py				      #
#		   italodovalle@gmail.com			      #
#######################################################################

import os
import sys
import argparse
from wes import ExomeFunctions as ef
import yaml
  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MuTect Analysis Pipeline')
    parser.add_argument('-n', required=True, dest='normal',
                        action='store',help='Normal bam file')
    parser.add_argument('-t', required=True, dest='tumor',
                        action='store',help='Tumor bam file')
    parser.add_argument('-s', required=True, dest='name',action='store',
                        help='Sample name')
    parser.add_argument('-c', required=True, dest='code', action='store',
                        help='Normal sample name')
    parser.add_argument('-f', required=True, dest='config', action='store',
                        help='config file yaml')
    parser.add_argument('-d', dest='outdir', action='store',
                        help='output directory', default= '.')
    if len(sys.argv) <= 3: 
        parser.print_help() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    

    ##--------------------------------
    #Input file fastq already trimmed
    normal = args.normal
    tumor = args.tumor
    name = args.name
    config = args.config
    outdir = os.path.abspath(args.outdir)
    code = args.code
    
    with open(config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
    
    outfile = '%s/%s_mutect.vcf'%(outdir,name)
    coverage_out = '%s/%s_coverage.wig'%(outdir,name)
    vcf = ef.run_mutect (cfg['softwares']['mutect'], 
                         cfg['ref-files']['hg'],
                         cfg['ref-files']['cosmic'],
                         cfg['ref-files']['dbsnp'],
                         cfg['sample-details']['target'],
                         normal, tumor, outfile,
                         coverage_out,outdir=outdir,err_log=True)
    mutect_dir = outdir+'/mutec_ann'
    os.mkdir(mutect_dir)
    sample_order = ef.get_order_vcf (code, vcf)
    ann = ef.Annotation(vcf,cfg['softwares']['convert2annovar'],
                        cfg['softwares']['annovar'],
                        cfg['ref-files']['humandb'],
                        cfg['ref-files']['build_ver'],
                        cfg['ref-files']['dbsnp_ver'],
                        cfg['ref-files']['kg_ver'],
                        cfg['ref-files']['mitochondrial_ver'],
                        fmt = 'vcf4old',outdir=mutect_dir,name=name,
                        err_log=True,maf=0.05,mutect=True)
    ef.MakeFinalFile (ann[0],ann[1],ann[2],
                      outfile='%s/%s.tsv'%(mutect_dir,name),
                      mutect=True,sample_order=sample_order)
