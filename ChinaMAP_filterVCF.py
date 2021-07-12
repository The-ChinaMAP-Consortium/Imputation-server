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
    sumFqList = []
    for sample in samples:
        try:
            GT = list(rec.samples.get(sample)['GT'])
            GT = [str(GT[0]), str(GT[1])]
            if "None" in GT:
                continue
            else:
                GT = GT[0] + GT[1]
                sumFqList .append(GT)
        except:
            continue
    if len(set(sumFqList)) == 1:
        return False
    else:
        return True

def out_site(rec, panelVcfFile):
    chrName = rec.chrom
    pos = rec.pos
    ref = rec.ref
    alts = list(rec.alts)
    if len(alts) > 1:
        return False
    getAlts = []
    for rec2 in panelVcfFile.fetch(chrName, pos-1, pos):   
        getPos =  rec2.pos
        getRef = rec2.ref
        getAlts += list(rec2.alts)
    try:
        if  (pos == getPos) and (ref == getRef) and (alts[0] in getAlts):
            return True
        else:
            return False
    except:
        return False
    
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
        if out_site(rec, panelVcfFile):
            pass
        else:
            outPanelNumber += 1
            outPanelFile.write(rec)
            continue
        outputFile.write(rec)
    return [monoNumber, outPanelNumber]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter the input files for monomorphic sites and none ChinaMAP imputation server sites.')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', type=str, default = None, required=True, nargs='+', help="Input .vcf.gz files.")
    required.add_argument('-d', '--outputDir', type=str, default = None, required=True, help="Output dictionary.")
    optional.add_argument('-r', '--reference', type=str, default = None, help="ChinaMAP reference panel .vcf.gz file with index file.")
    optional.add_argument('-e', '--exclude_monomorphic',  action='store_true', help="Exclude monomorphic sites.")
    args = parser.parse_args()
    inFileList = args.input
    outDir = args.outputDir
    
    LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename='ChinaMAP_filterVCF.log', level=logging.DEBUG, format=LOG_FORMAT)

    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    for fileName in inFileList:
        if fileName.endswith(".vcf.gz"):
            try:
                tmpFile = VariantFile(fileName, 'r')
                tmpFile.close()            
            except:
                logging.error(fileName + ": Unsupported file types gunzip.")
                continue
            else:
                prefix = fileName.strip('.vcf.gz').split('/')[-1]
        else:
            logging.error(vcfFile.filename.decode() + ": Not in .vcf.gz format.")
            continue

        if not (args.reference or args.exclude_monomorphic):
            logging.error('No action requested, add -r or -e')
        
        outName = outDir+"/"+prefix+"_filtered.vcf.gz"
        vcfFile = VariantFile(fileName, 'r')

        outFile = VariantFile(outName, 'w', header = vcfFile.header)

        if (args.reference and args.exclude_monomorphic):

            logging.info(vcfFile.filename.decode() + ": Start filtering monomorphic sites and none ChinaMAP reference panel sites.")

            referencePanel = args.reference
            panelVcfFile = VariantFile(referencePanel, 'r')
            outPanelName = outDir+"/"+prefix+"_exclude.vcf.gz"
            monoName = outDir+"/"+prefix+"_monomorphic.vcf.gz"

            outPanelFile = VariantFile(outPanelName, 'w', header = vcfFile.header)
            monoFile   = VariantFile(monoName, 'w', header = vcfFile.header)
            
            numbers = vcf_filter(vcfFile, panelVcfFile, monoFile, outPanelFile, outFile)
            
            vcfFile.close()
            outFile.close()
            outPanelFile.close()
            monoFile.close()

            logging.info(vcfFile.filename.decode() + ": Excluded {monoNumber} monomorphic sites and {outPanelNumber} none ChinaMAP reference panel sites.".format(monoNumber = numbers[0], outPanelNumber = numbers[1]))
        
        elif args.reference:

            logging.info(vcfFile.filename.decode() + ": Start filtering none ChinaMAP reference panel sites.")

            referencePanel = args.reference
            panelVcfFile = VariantFile(referencePanel, 'r')
            outPanelName = outDir+"/"+prefix+"_exclude.vcf.gz"

            outPanelFile = VariantFile(outPanelName, 'w', header = vcfFile.header)

            outPanelNumber = 0

            for rec in vcfFile.fetch():
                if out_site(rec, panelVcfFile):
                    outFile.write(rec)
                else:
                    outPanelNumber += 1
                    outPanelFile.write(rec)

            vcfFile.close()
            outFile.close()
            outPanelFile.close()

            logging.info(vcfFile.filename.decode() + ": Excluded {outPanelNumber} none ChinaMAP reference panel sites.".format(outPanelNumber = outPanelNumber))

        else: 

            logging.info(vcfFile.filename.decode() + ": Start filtering monomorphic sites.")

            monoName = outDir+"/"+prefix+"_monomorphic.vcf.gz"
            
            monoFile   = VariantFile(monoName, 'w', header = vcfFile.header)
            
            monoNumber = 0
            samples = list(vcfFile.header.samples)
            for rec in vcfFile.fetch():
                if mono_site(rec, samples):
                    outFile.write(rec)
                else:
                    monoNumber += 1
                    monoFile.write(rec)
                
            vcfFile.close()
            outFile.close()
            monoFile.close()

            logging.info(vcfFile.filename.decode() + ": Excluded {monoNumber} monomorphic sites.".format(monoNumber = monoNumber))