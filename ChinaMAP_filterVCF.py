#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   ChinaMAP_filterVCF.py
@Time    :   2021/07/09 13:29:14
@Author  :   Siyu Wang
@Version :   1.0
@Contact :   wangsiyu.13@sjtu.edu.cn
@Desc    :   None
'''

import argparse
import os
from pysam import VariantFile
import logging

def mono_site(rec, samples):
    sumFqList = [sum(rec.samples.get(sample)['GT']) for sample in samples]
    if len(set(sumFqList)) == 1:
        return True
    else:
        return False

def out_site(rec, panelVcfFile):
    chrName = rec.chrom
    pos = rec.pos
    for rec in panelVcfFile.fetch(chrName, pos, pos+1):
        panelPos = rec.pos
    
    if panelPos == pos:
        return False
    else:
        return True
    
def vcf_filter(vcfFile, panelVcfFile, monoFile, outPanelFile, outputFile):
    monoNumber = 0
    outPanelNumber = 0
    samples = list(vcfFile.header.samples)
    for rec in vcfFile.fetch():
        if mono_site(rec, samples):
            pass
        else:
            monoNumber += 1
            monoFile.write(rec)
            continue
        if out_site(rec):
            pass
        else:
            outPanelNumber += 1
            outPanelFile.write(rec)
            continue
        outputFile.write(rec)
    return [monoNumber, outPanelNumber]

def write_header(tplVcfFile, otherFile):
    for rec in vcfFile.header.records:
        if rec.key in ['fileformat', 'FILTER', 'ALT']:
            otherFile.header.add_record(rec)
    otherFile.header.add_record(vcfFile.header.formats.get('GT').record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check input file for ChinaMAP imputation server.')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', type=str, default = None, required=True, help="Input .vcf.gz files path list.")
    required.add_argument('-d', '--outputDir', type=str, default = None, required=True, help="Output dictionary.")
    optional.add_argument('-r', '--reference', type=str, default = None, help="ChinaMAP reference panel .vcf.gz file.")
    optional.add_argument('-e', '--exlude_monomorphic',  action='store_true', help="Exclude monomorphic sites.")
    optional.add_argument('-m', '--minimize', action='store_true', help="Delete non-essential information from .vcf.gz file.")
    args = parser.parse_args()
    inFileList = args.input
    outDir = args.outputDir
    
    LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename='ChinaMAP_filterVCF.log', level=logging.DEBUG, format=LOG_FORMAT)

    
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    fileList = []
    with open(inFileList, "r") as f:
        for line in f.readlines():
            fileList.append(line.strip('\n'))
    
    for fileName in fileList:
        if fileName.endswith(".vcf.gz"):
            prefix = fileName.strip('.vcf.gz').split('/')[-1]
        else:
            print()
            continue

        if not (args.reference or args.exclude_monomorphic):
            logging.error('No action requested, add -r or -e')
        
        outName = "outDir"+"/"+"prefix"+"_filtered.vcf.gz"
        vcfFile = VariantFile(fileName, 'r')
        outFile = VariantFile(outName, 'w')

        if args.minimize:
            for rec in vcfFile.fetch():
                rec.info.clear
                samples = list(vcfFile.header.samples)
                for sample in samples:
                    GT = rec.samples.get(sample)['GT']
                    rec.samples.clear()
                    rec.samples.get(sample)['GT'] = GT        

        if (args.reference and args.exclude_monomorphic):

            logging.info(vcfFile.filename + ": Start filtering monomorphic sites and none ChinaMAP reference panel sites.")

            referencePanel = args.reference
            panelVcfFile = VariantFile(referencePanel, 'r')
            outPanelName = "outDir"+"/"+"prefix"+"_exclude.vcf.gz"
            monoName = "outDir"+"/"+"prefix"+"_monomorphic.vcf.gz"

            if args.minimize:
                outPanelFile = VariantFile(outPanelName, 'w')
                monoFile   = VariantFile(monoName, 'w')
                write_header(vcfFile, outPanelFile)
                write_header(vcfFile, monoFile)
                write_header(vcfFile, outFile)
            else:
                outPanelFile = VariantFile(outPanelName, 'w', header = vcfFile.header)
                monoFile   = VariantFile(monoName, 'w', header = vcfFile.header)
            
            numbers = vcf_filter(vcfFile, panelVcfFile, monoFile, outPanelFile, outFile)
            
            logging.info(vcfFile.filename + ": Excluded {monoNumber} monomorphic sites and {outPanelNumber} none ChinaMAP reference panel sites.".format(monoNumber = numbers[0], outPanelNumber = numbers[1]))
        
        elif args.reference:

            logging.info(vcfFile.filename + ": Start filtering none ChinaMAP reference panel sites.")

            referencePanel = args.reference
            panelVcfFile = VariantFile(referencePanel, 'r')
            outPanelName = "outDir"+"/"+"prefix"+"_exclude.vcf.gz"

            if args.minimize:
                outPanelFile = VariantFile(outPanelName, 'w')
                write_header(vcfFile, outPanelFile)
                write_header(vcfFile, outFile)
            else:
                outPanelFile = VariantFile(outPanelName, 'w', header = vcfFile.header)

            outPanelNumber = 0

            for rec in vcfFile.fetch():
                if out_site(rec, panelVcfFile):
                    outFile.write(rec)
                else:
                    outPanelNumber += 1
                    outPanelFile.write(rec)

            logging.info(vcfFile.filename + ": Excluded {outPanelNumber} none ChinaMAP reference panel sites.".format(outPanelNumber = outPanelNumber))

        else: 

            logging.info(vcfFile.filename + ": Start filtering monomorphic sites.")

            monoName = "outDir"+"/"+"prefix"+"_monomorphic.vcf.gz"
            
            if args.minimize:
                monoFile   = VariantFile(monoName, 'w')
                write_header(vcfFile, monoFile)
                write_header(vcfFile, outFile)
            else:
                monoFile   = VariantFile(monoName, 'w', header = vcfFile.header)
            
            monoNumber = 0
            samples = list(vcfFile.header.samples)
            for rec in vcfFile.fetch():
                if mono_site(rec, samples):
                    outFile.write(rec)
                else:
                    monoNumber += 1
                    monoFile.write(rec)

            logging.info(vcfFile.filename + ": Excluded {monoNumber} monomorphic sites.".format(monoNumber = monoNumber))