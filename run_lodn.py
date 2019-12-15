#! /usr/bin/env python

#######################################################################
#		     Italo Faria do Valle                             #
#								      #
#		        run_lodn.py				      #
#		   italodovalle@gmail.com			      #
#######################################################################

import os
from wes import lodn
import argparse
import sys

if __name__ == '__main__':
    description = 'Run LODN filtering'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', required=True, dest='infile',action='store',
                        help='infile')
    parser.add_argument('-f', required=True, dest='fmt',action='store',
                        help='table or vcf', default='table')
    parser.add_argument('-b', required=True, dest='bam', action='store',
                        help='normal bam file')
    parser.add_argument('-o', required=True, dest='outfile', action='store',
                        help='outfile')
    parser.add_argument('-min_n_cov', dest='min_n_cov', action='store',
                        help='min normal coverage', default=8)
    parser.add_argument('-min_t_cov', dest='min_t_cov', action='store',
                        help='min tumor coverage', default=14)
    
    if len(sys.argv) <= 3: 
        parser.print_help() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    ##--------------------------------
    infile = args.infile
    fmt = args.fmt
    bam = args.bam
    outfile = args.outfile
    min_n_cov = int(args.min_n_cov)
    min_t_cov = int(args.min_t_cov)
    
        
    if fmt == 'table':
        lodn.lodn_calculator_table(infile, '\t', bam, outfile=outfile,
                                   min_n_cov = min_n_cov, 
                                   min_t_cov = min_t_cov,
                                   lodn_known = 5.5, lodn_novel=2.2)
    elif fmt == 'vcf':
        lodn.lodn_calculator_vcf(infile, bam, outfile=outfile,
                                 min_n_cov = min_n_cov,
                                 min_t_cov = min_t_cov,
                                 lodn_known = 5.5, lodn_novel=2.2)