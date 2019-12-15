#! /usr/bin/env python

#######################################################################
#		     Italo Faria do Valle                             #
#								      #
#		      ExomeFunctions.py				      #
#		   italodovalle@gmail.com			      #
#######################################################################

import os
import sys
import time
import yaml
import subprocess
import re
import pandas as pd
from pandas import DataFrame
from collections import defaultdict

def verify (file,func_name):
    if (os.path.getsize(("%s"%file))) == 0:
        msg = "\n\nERROR: File %s empty! Check %s step!\n\n"
        sys.exit((msg%(file,func_name)))


def create_tree_dir(local='.'):
    local = os.path.abspath(local)
    """Create directory tree for organize results"""
    analysis = ['/fastq','/Coverage', '/Data-cleanup','/Variant-Discovery',
                '/Results']
    data_cleanup = ['/Aligned', '/BQRS', '/IndelRealigned', '/Logs',
                    '/Marked', '/Sorted', '/Trimmed']
    os.mkdir(local+'/Analysis')
    local = local+'/Analysis'
    for i in analysis:
        os.mkdir(local+i)
    for i in data_cleanup:
        os.mkdir(local+'/Data-cleanup'+i)
    os.mkdir(local+'/Variant-Discovery/Logs')
    
def Mapping(fastq,bwa,ref,outdir='.',name='sample',PL='platform',
            LB='library',cores=2,err_log = True):
    """Mapping to the human reference"""
    outdir = os.path.abspath(outdir)
    
    ## Breaking command line
    inputs = len(fastq)*'%s '%(tuple(fastq))
    output = '%s/%s.sam'%(outdir, name)
    RG = '\'@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\\tPU:PE\''%(name,name,PL,LB)
    basic = '%s mem -t %d -M -R %s %s '%(bwa, cores, RG, ref)
    if err_log:
        log = '2>%s/err.bwa_run'%outdir
    else:
        log = ''
    os.system(basic+inputs+'> '+output +' '+ log)
    verify(output,Mapping.__name__)
    
    return (output)

def Sort(sam,picard, tmp_dir='/tmp',outdir='.',name='sample',err_log = True):
    """Sorting alignment """
    outdir = os.path.abspath(outdir)
    sort = 'java -Djava.io.tmpdir=%s -jar  %s SortSam '%(tmp_dir,picard)
    output = '%s/%s_sorted.bam'%(outdir, name)
    parameters = 'INPUT=%s OUTPUT=%s SORT_ORDER=coordinate '%(sam, output)
    if err_log:
        log = '2>%s/err.SortSam'%outdir
    else:
        log = ''
    os.system(sort+parameters+'TMP_DIR=%s'%tmp_dir+log)
    verify(output, Sort.__name__)
    
    return (output)

def MarkDuplicates(bam,picard, tmp_dir='/tmp',outdir='.',name='sample',
                   err_log=True,
                   MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=None):
    
    outdir = os.path.abspath(outdir)
    java = 'java -Djava.io.tmpdir=%s -jar'%tmp_dir
    mark = '%s MarkDuplicates'%(picard)
    output = '%s/%s_dedup.bam'%(outdir, name)
    metrics = 'METRICS_FILE=%s/metrics.txt'%outdir
    parameters = 'INPUT=%s OUTPUT=%s %s'%(bam, output, metrics)
    if err_log:
        log='2>%s/err.MarkDuplicates'%outdir
    else:
        log = ''
        
    ## when the picard raises an error of not finding the tmp files
    if MAX_FILE_HANDLES_FOR_READ_ENDS_MAP:
        subprocess.call(' '.join([java, mark,parameters,
                                  'TMP_DIR=%s'%tmp_dir,
                                  'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=',
                                  str(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP),
                                  log]),
                        shell=True)
    else:
        subprocess.call(' '.join([java, mark,parameters,
                                  'TMP_DIR=%s'%tmp_dir,
                                  log]),
                        shell=True)
    verify(output, MarkDuplicates.__name__)
    
    return (output)

def BuildBamIndex (bam, picard, tmp_dir='/tmp', outdir='.', name='sample', 
                   err_log=True):
    outdir=os.path.abspath(outdir)
    java = 'java -Djava.io.tmpdir=%s -jar '%(tmp_dir)
    build = '%s BuildBamIndex'%(picard)
    parameters = 'INPUT=%s TMP_DIR=%s'%(bam,tmp_dir)
    if err_log:
        log='2>%s/err.BuildBamIndex'%outdir
    else:
        log = ''
    command= java+' '+build+' '+parameters+' '+log
    subprocess.call(command, shell=True)

def GetUnmapped (bam,picard,name='unmapped', outdir='.',err_log=True):
    """Get unmapped reads from bam file"""
    ## Pay Attention to the Picard Version
    outdir=os.path.abspath(outdir)
    f1 = '%s/%s_1.fastq'%(outdir,name)
    f2 = '%s/%s_2.fastq'%(outdir,name)
    fu = '%s/%s_unpaired.fastq'%(outdir,name)
    samtofastq='java -jar %s SamToFastq'%picard
    files='I=%s F=%s F2=%s FU=%s'%(bam,f1,f2,fu)
    command = samtofastq+files#+'VALIDATION_STRINGENCY=LENIENT'
    if err_log:
    	log='2>%s/err.SamToFastq'%outdir
    else:
        log=''
    subprocess.call(' '.join([samtofastq,files,log]), shell=True)
    
    return ([f1, f2, fu])
    
def IndelRealigner (bam, gatk, ref, target, indels_ref,outdir='.',
                    name='sample',err_log=True,tmp_dir='/tmp', L=None, 
                    ip=None):
    """Local Realigment around indels"""
    outdir = os.path.abspath(outdir)
    if err_log:
        log='2>%s/err.IndelRealigner'%outdir
    else:
        log = ''
    
    java = 'java -Djava.io.tmpdir=%s -jar'%(tmp_dir)
    cmd = '%s -T IndelRealigner'%(gatk)
    output = '%s/%s_realigned.bam'%(outdir,name)
    
    command = ' '.join([java,cmd, '-R', ref, '-I', bam, 
                              '-targetIntervals',target,'-known',indels_ref,
                              '-o', output])
    
    if L:
        command = ' '.join([command, '-L', L])
    if ip:
        command = ' '.join([command, '-ip', str(ip)])
    
    command = ' '.join([command, log])
    
    subprocess.call(command, shell=True)
    verify(output, IndelRealigner.__name__)
    
    return(output)


def BQSR(bam, gatk, ref,dbsnp, indels_ref,name='sample',cores=2, outdir='.',
         err_log=True,tmp_dir='/tmp',L=None,ip=None):
    """Base Quality Score Recalibration"""
    cores = str(cores)
    outdir = os.path.abspath(outdir)
    if err_log:
        log = ['2>%s/err.BaseRecalibrator'%outdir,
                '2>%s/err.BaseRecalibrator_2round'%outdir,
                '2>%s/err.AnalyseCovariates'%outdir,
                '2>%s/err.PrintReads'%outdir]
    else:
        log = ['','','','']
    
    
    java = 'java -Djava.io.tmpdir=%s -jar %s'%(tmp_dir,gatk)
    outtable1 = '%s/%s_recal_data.table'%(outdir,name)
    
    command = ' '.join([java, '-T BaseRecalibrator','-R', ref,'-I',
                              bam,'-knownSites',dbsnp,'-knownSites',
                              indels_ref,'-nct',cores,'-o',outtable1])
    if L:
        command = ' '.join([command, '-L', L])
    if ip:
        command = ' '.join([command, '-ip', str(ip)])

    command = ' '.join([command, log[0]])
    
    subprocess.call(command, shell=True)
    
    outtable2 = '%s/%s_post_recal_data.table'%(outdir,name)
    subprocess.call(' '.join([java, '-T BaseRecalibrator','-R', ref,'-I',
                              bam,'-knownSites',dbsnp,'-knownSites',
                              indels_ref,'-BQSR', outtable1, '-nct',cores,
                              '-o',outtable2,log[1]]), shell=True)
    
    plots =  '%s/%s_recalibrationPlots.pdf'%(outdir,name)

    subprocess.call(' '.join([java, '-T AnalyzeCovariates', '-R', ref,
                              '-before', outtable1, '-after', outtable2,
                              '-plots', plots,log[2]]), shell=True)
    
    recal = '%s/%s_recal.bam'%(outdir,name)
    subprocess.call(' '.join([java, '-T PrintReads', '-R', ref, '-I', bam,
                             '-BQSR', outtable1, '-nct', cores, '-o', 
                             recal, log[3]]), shell=True)
   
    verify(recal, BQSR.__name__)
    return(recal)

 
def HaplotypeCaller(bam,gatk,ref,outdir='.',name='sample',err_log=True,
                    genotyping_mode = 'DISCOVERY',stand_emit_conf=10,
                    stand_call_conf=30):  
    """GATK Haplotype Caller"""
    outdir = os.path.abspath(outdir)
    hc = 'java -jar %s -T HaplotypeCaller'%gatk
    output = '%s/%s_raw.vcf'%(outdir,name)
    if err_log:
        log='2>%s/err.HaplotypeCaller'%outdir
    else:
        log = ''
    
    cmd = ' '.join([hc, '-R', ref, '-I', bam, '--genotyping_mode', 
                    genotyping_mode, '-stand_emit_conf', str(stand_emit_conf),
                    '-stand_call_conf', str(stand_call_conf),'-o', output,
                    log])
    
    subprocess.call(cmd, shell=True)
    
    verify(output,HaplotypeCaller.__name__)
    return(output)
    
def VQSR_SNP (vcf,gatk,ref,res_hapmap,res_omni, omni, res_dbsnp, dbsnp,
              snp_pars,ts_filter_level,outdir='.',name='sample',
              err_log=True):
    """GATK Variant Quality Score Recalibration SNPs"""
    outdir = os.path.abspath(outdir)
    
    ## Build the SNP recalibration model
    var_recal = 'java -jar %s -V VariantRecalibrator -R %s'%(gatk,ref)
    inp = '-input %s'%(vcf)
    resource = 3*'-resource:%s %s '%(res_hapmap, hapmap,res_omni, omni,
                                     res_dbsnp, dbsnp)
    recalFile = '-recalFile %s/%s_recalibrate_SNP.recal'%(outdir,name)
    tranchesFile = '-tranchesFiles %s/%s_recalibrate_SNP.tranches'%(name,
                                                                    outdir)
    rscriptFile = '-rscriptFile %s/%s_recalibrate_SNP_plots.R'%(name,outdir)
    
    if err_log:
        log='2>%s/err.VariantRecalibrator_SNP'%outdir
    else:
        log = ''
    
    subprocess.call(' '.join([var_recal,inp,resource,snp_pars,recal_file,
                              tranchesFile,rscriptFile,log]))

    #2) Apply the desired level of recalibration to the SNPs in the call set
    applyrecal = 'java -jar %s -T ApplyRecalibration -R %s'%(gatk,ref)
    par = '-mode SNP --ts_filter_level %f'%float(ts_filter_level)
    out = '-o %s/%s_recalibrated_snps.vcf'%(outdir,name)
    
    if err_log:
        log='2>%s/err.ApplyRecalibration_SNP'%outdir
    else:
        log = ''

    subprocess.call(' '.join([applyrecal,inp,par,recalFile,tranchesFile,out,
                              log]), shell=True)
    verify(out,VQSR_SNP.__name__)
    
def VQSR_INDEL (vcf,gatk,ref,res_mills,indels_ref,indels_pars,
                ts_filter_level,outdir='.',name='sample',
                err_log=True):
    """GATK Variant Quality Score Recalibration INDELs"""
    outdir = os.path.abspath(outdir)
    
    ## Build the INDEL recalibration model
    var_recal = 'java -jar %s -V VariantRecalibrator -R %s'%(gatk,ref)
    inp = '-input %s'%(vcf)
    resource = '-resource:%s %s '%(res_mills,indels_ref)
    recalFile = '-recalFile %s/%s_recalibrate_INDEL.recal'%(outdir,name)
    tranchesFile = '-tranchesFiles %s/%s_recalibrate_INDEL.tranches'%(name,
                                                                    outdir)
    rscriptFile = '-rscriptFile %s/%s_recalibrate_INDEL_plots.R'%(name,outdir)
    
    if err_log:
        log='2>%s/err.VariantRecalibrator_INDEL'%outdir
    else:
        log = ''
        
    subprocess.call(' '.join([var_recal, inp, resource, indels_pars,
                              recalFile,tranchesFile,rscriptFile,log]))

    ## Apply the desired level of recalibration to the Indels in the call set
    applyrecal = 'java -jar %s -T ApplyRecalibration -R %s'%(gatk,ref)
    par = '-mode INDEL --ts_filter_level %f'%float(ts_filter_level)
    out = '%s/%s_recalibrated_indels.vcf'%(outdir,name)
    
    if err_log:
        log='2>%s/err.ApplyRecalibration_INDEL'%outdir
    else:
        log = ''

    subprocess.call(' '.join([applyrecal,inp,par,recalFile,tranchesFile,'-o',
                              out,log]), shell=True)
    verify(out,VQSR_INDEL.__name__)

def HardFilter (vcf,gatk,ref,filter_exp_snps, filter_exp_indels,outdir='.',
                name='sample',err_log=True):
    """GATK Hard Filtering Variants"""
    selectVars = 'java -jar %s -T SelectVariants -R %s -V'%(gatk, ref)
    
    ## Extract the snps from the call set
    selectType = '-selectType SNP'
    raw_snps = '%s/%s_snps.vcf'%(outdir,name)
    if err_log:
        log='2>%s/err.SelectSnpsHF'%outdir
    else:
        log = ''
    subprocess.call(' '.join([selectVars,vcf,selectType,'-o',raw_snps,log]),
                    shell=True)
    
    ## Apply the filter to the SNP call set
    var_recal = 'java -jar %s -T VariantFiltration'%(gatk)
    filt_snps = '%s/%s_hard_filtered_snps.vcf'%(outdir,name)
    if err_log:
        log='2>%s/err.applyHFSNPs'%outdir
    else:
        log = ''
    subprocess.call(' '.join([var_recal,'-R', ref, '-V', raw_snps, 
                              filter_exp_snps,'-o',filt_snps,
                              log]), shell=True)
    
    ## Extract the Indels from the call set
    selectType = '-selectType INDEL'
    raw_indels = '%s/%s_indels.vcf'%(outdir,name)
    if err_log:
        log='2>%s/err.SelectIndelsHF'%outdir
    else:
        log = ''
    subprocess.call(' '.join([selectVars,vcf,selectType,'-o',raw_indels,log]),
                    shell=True)
    
    ## Apply the filter to the Indel call set
    if err_log:
        log='2>%s/err.applyHFIndels'%outdir
    else:
        log = ''
    filt_indels = '%s/%s_hard_filtered_indels.vcf'%(outdir,name)
    subprocess.call(' '.join([var_recal,'-R', ref, '-V', raw_indels,
                              filter_exp_indels,'-o',
                              filt_indels,log]),shell=True)
    if err_log:
        log='2>%s/err.CombineVariants'%outdir
    else:
        log = ''
    
    ## Combine variants
    combine = 'java -jar %s -T CombineVariants -R %s'%(gatk,ref)
    out = '%s/%s_filtered_variants.vcf'%(outdir,name)
    subprocess.call(' '.join([combine,'--variant:snps',filt_snps,
                              '--variant:indels',filt_indels,'-o',out,
                              '-genotypeMergeOptions PRIORITIZE',
                              '-priority snps,indels' ,log]),
                    shell=True)
    
    verify(out,HardFilter.__name__)

    return (out)

def annovar_filter(annovar,infile, build_ver,dbtype, humandb,outdir='.',
                   err_log=True, pars=''):
    
    if err_log:
        log='2>%s/err.annovar_%s'%(outdir,dbtype)
    else:
        log = ''
    
    subprocess.call(' '.join([annovar, '-filter', infile,'-buildver',
                              build_ver,'-dbtype', dbtype, humandb,
                              pars,log]),shell=True)
    
def process_annovar_out(dbtype, ann, rmdup,mutect=False):
    if mutect:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$17,$18,$19}\''%dbtype
    else:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$19,$20}\''%dbtype
    subprocess.call(' '.join([command, ann,'>', rmdup]),shell=True)

def Annotation(vcf,convert2annovar,annovar,humandb,build_ver, dbsnp_ver,
               kg_ver,mitochondrial_ver,fmt = 'vcf4',outdir='.',name='sample',
               err_log=True,maf=0.05,mutect=False):
    
    if mutect:
        fmt = 'vcf4old'
        pars = '--filter pass'
    else:
        pars = ''
    
    # convert to annovar format
    outdir = os.path.abspath(outdir)
    outfile = '%s/%s.annovar'%(outdir,name)
    
    if err_log:
        log='2>%s/err.convert2annovar'%outdir
    else:
        log = ''
    subprocess.call(' '.join([convert2annovar, '-format',fmt, vcf, '-outfile',
                              outfile, '-includeinfo', '-withzyg',
                              '-comment', pars,log]),shell=True)
    # annotate
    #dbSNP138 and 1000g annotation
    annovar_filter(annovar, outfile,build_ver, dbsnp_ver,humandb,outdir)
    dbsnp_ann = '%s.%s_%s_dropped'%(outfile,build_ver,dbsnp_ver)
    dbsnp_rmdup = '%s/%s_rmdup.dbsnp'%(outdir,name)
    process_annovar_out(dbsnp_ver, dbsnp_ann, dbsnp_rmdup,mutect)
    dbsnp_filt = '%s.%s_%s_filtered'%(outfile,build_ver,dbsnp_ver)
    
    annovar_filter(annovar, dbsnp_filt, build_ver,kg_ver,humandb,outdir,
                   pars='-maf 0.05 -reverse')
    suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,
                                                    kg_ver[-3:].upper())
    
    kg_ann = [x for x in os.listdir(outdir) if re.findall(suffix,x)][0]
    kg_ann = outdir+'/'+kg_ann
    kg_rmdup = '%s/%s_rmdup.1000g'%(outdir,name)
    process_annovar_out(kg_ver, kg_ann, kg_rmdup,mutect)
    
    known_file = '%s/%s_rmdup.known'%(outdir,name)
    novel_file = '%s/%s_rmdup.novel'%(outdir,name)
    
    os.system('cat %s %s > %s'%(dbsnp_rmdup,kg_rmdup,known_file))
    
    #Turning your life more easy =) -> one changing the name of novel variants file
    os.system('mv %s %s'%(kg_ann,novel_file))
    
    ## Gene-based annotation
    for ann_file in [known_file, novel_file]:
        subprocess.call(' '.join([annovar,'-geneanno',ann_file,'-buildver',
                              build_ver,humandb]), shell=True)
    
    #Mitochondrial Annotation
    
    subprocess.call(' '.join([annovar, '-buildver', mitochondrial_ver,
                              '-dbtype ensGene', outfile, humandb,
                              log]), shell=True)
    if mutect:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$19,$20,$21}\''
        command = command%mitochondrial_ver
    else:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$21,$22,$23}\''
        command = command%mitochondrial_ver
    mit_ann = '%s.exonic_variant_function'%outfile
    mit_rmdup = '%s/%s.mit'%(outdir,name)
    subprocess.call(' '.join([command, mit_ann, '>',mit_rmdup]),shell=True)
    
    ann_files = ['%s.exonic_variant_function'%known_file,
                 '%s.exonic_variant_function'%novel_file,
                 mit_rmdup]
    
    return(ann_files)

def parse_known (fi_known, mutect=False,sample_order=['n','t']):
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    
    cmd = "cat %s | tr [:blank:] '\t'"%fi_known
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE,shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    known = output.split('\n')
    if known[-1] == '':
        known = known[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info_%s'%sample_order[0],
                'info_%s'%sample_order[1]]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info']
    for i in range(len(known)-1):
        vec = known[i].split('\t')
        for k in range(len(vec)):
            dic[i][cols[k]] = vec[k]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        gt = [x.split(':')[0] for x in dt['info_t']]
        ## ugly solution
        ## bug: different genotype coding found in vcf. The code was set for
        ## genotypes coded as in the genotype dictionary
        try:
            dt['genotype'] = [genotype[x] for x in gt]
        except:
            dt['genotype'] = float('nan')

        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]               
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['known.flag'] = 1
    return(dt)

def parse_novel (infile,mutect=False,sample_order=['n','t']):
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    cmd = "cat %s | tr [:blank:] '\t'"%infile
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    data = output.split('\n')
    if data[-1] == '':
        data = data[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info_%s'%sample_order[0],
                'info_%s'%sample_order[1]]
        for i in range(len(data)):
            vec = data[i].split('\t')
            for k in range(len(vec[:11])):
                dic[i][cols[k]] = vec[k]
            dic[i]['format'] = vec[-3]
            dic[i]['info_n'] = vec[-2]
            dic[i]['info_t'] = vec[-1]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info']
        for i in range(len(data)-1):
            vec = data[i].split('\t')
            for k in range(len(vec[:11])):
                dic[i][cols[k]] = vec[k]
            dic[i]['format'] = vec[-2]
            dic[i]['info'] = vec[-1]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        gt = [x.split(':')[0] for x in dt['info_t']]
        ## ugly solution
        ## bug: different genotype coding found in vcf. The code was set for
        ## genotypes coded as in the genotype dictionary
        try:
            dt['genotype'] = [genotype[x] for x in gt]
        except:
            dt['genotype'] = float('nan')
        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]        
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['name_var'] = 'novel'
    dt['known.flag'] = 0
    return(dt)
    
   
def parse_mit (fi_known,mutect=False,sample_order=['n','t']):  
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    cmd = "cat %s | tr [:blank:] '\t'"%fi_known
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE,shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    data = output.split('\n')
    ## check if last line is empty
    if data[-1] == '':
        data = data[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','database','type',
                'format','info_%s'%sample_order[0], 
                'info_%s'%sample_order[1]]
        for i in range(len(data)):
            vec = data[i].rstrip().split('\t')
            #ixs = [2,3,4,5,6,7,8,10,1,18,19,20]
            ixs = [0,1,2,3,4,5,6,7,8,9,10,11]
            for k in range(len(ixs)):
                dic[i][cols[k]] = vec[ixs[k]]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','database','type',
                'format','info']
        for i in range(len(data)):
            vec = data[i].rstrip().split('\t')
            for k in range(len(vec)):
                dic[i][cols[k]] = vec[k]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]        
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['genotype'] = 'Mit'
    dt['name_var'] = 'Mit'
    dt['known.flag'] = 0
    return(dt)

def MakeFinalFile (fi_known,fi_novel,fi_mit,outfile='sample.tsv',mutect=False,
                   sample_order=['n','t'],dbsnp_freq=True, dbsnpFreq = None,
                   dbsnpAllele=None):
    
    ## defining column names
    cols = []
    if mutect:
        cols = ['name_var','type','chr', 'pos.start', 'pos.end','ref.base',
                'alt.base','genotype','annotation','t_cov.ref','t_cov.alt',
                'n_cov.ref','n_cov.alt','known.flag']
    else:
        cols = ['name_var','type','chr', 'pos.start', 'pos.end','ref.base',
                'alt.base','genotype','annotation','cov.ref','cov.alt',
                'known.flag']
    
    test_size = 0
    for i in [fi_novel, fi_known, fi_mit]:
        test_size += os.path.getsize(i)
    if test_size == 0:
        final = DataFrame(columns=cols)
    else:
        dt = {}
        if os.path.getsize(fi_novel) > 0:
            dt['novel'] = parse_novel(fi_novel,mutect=mutect,
                                      sample_order=sample_order)
            dt['novel']['known.flag'] = 0
        if os.path.getsize(fi_known) > 0:
            dt['known'] = parse_known(fi_known,mutect=mutect,
                                      sample_order=sample_order)
            dt['known']['known.flag'] = 1
        if os.path.getsize(fi_mit) > 0:
            dt['mit'] = parse_mit(fi_mit,mutect=mutect,
                                  sample_order=sample_order)
            dt['mit']['known.flag'] = 0
        for i in dt.keys():
            if mutect:
                dt[i]['genotype'] = 'unknown'
            dt[i] = dt[i][cols]
        
        final = pd.concat(dt.values(), ignore_index=True)
        
        if dbsnp_freq:
            if dbsnpFreq != None and dbsnpAllele != None:
                ## prepare dict
                dbsnp = defaultdict(dict)
                for i in open(dbsnpFreq).readlines():
                    vec = i.rstrip().split('\t')
                    dbsnp[vec[0]]['allele_id'] = vec[1]
                    dbsnp[vec[0]]['freq'] = vec[3]
                    dbsnp[vec[0]]['allele'] = ''
                    dbsnp[vec[0]]['allele_rev'] = ''
                alleles = defaultdict(dict)
                for i in open(dbsnpAllele).readlines():
                    vec = i.rstrip().split('\t')
                    alleles[vec[0]]['allele'] = vec[1]
                    alleles[vec[0]]['allele_rev'] = vec[3]
                for rs in dbsnp.keys():
                    try:
                        id1 = dbsnp[rs]['allele_id']
                        dbsnp[rs]['allele'] = alleles[id1]['allele']
                        id2 = alleles[dbsnp[rs]['allele_id']]['allele_rev']
                        dbsnp[rs]['allele_rev'] = alleles[id2]['allele']
                    except:
                        next
                        
            ## check freq
            final['dbsnp.freq'] = float('NaN')
            for i in final.index:
                var = final['name_var'].loc[i]
                if var.startswith('rs'):
                    try:
                        final['dbsnp.freq'].loc[i] = dbsnp_dic[var[2:]]
                    except:
                        next
        
    final.to_csv(outfile, sep='\t',header=True, index=None)

def run_mutect (mutect, ref, cosmic, dbsnp, intervals,normal, tumor, outfile,
                coverage_out,outdir='.',err_log=True):
    outdir = os.path.abspath(outdir)
    if err_log:
        log='1>%s/err.mutect'%outdir
    else:
        log = ''
    
    subprocess.call(' '.join([mutect, '--analysis_type MuTect',
                              '--reference_sequence',ref,'--cosmic',
                              cosmic, '--intervals',intervals,
                              '--input_file:normal', normal,
                              '--input_file:tumor', tumor,
                              '--vcf',outfile,'--coverage_file',
                              coverage_out,log]), shell=True)
    verify(outfile, run_mutect.__name__)
    return(outfile)

def get_order_vcf (normal_name, infile):
    
    for line in open(infile,'r').readlines():
        if line.startswith('#CHROM'):
            order = line.rstrip().split('\t')[-2:]
            if order.index(normal_name) == 1:
                return(['t','n'])
            if order.index(normal_name) == 0:
                return(['n','t'])