#!/usr/bin/env python3
# -*-coding:utf-8 -*-
'''
@File    :   ChinaMAP_checkVCF.py
@Time    :   2021/07/07 10:32:14
@Author  :   Siyu Wang
@Version :   1.0
@Contact :   wangsiyu.13@sjtu.edu.cn
@Desc    :   None
'''
import argparse
import os
from pysam import VariantFile
import logging

def get_FileSize(filePath):
    fsize = os.path.getsize(filePath)
    fsize = fsize/float(1024*1024)
    if fsize > 100:
        logging.info(filePath + ": File size over 100M.")
        return False
    else:
        logging.info(filePath + ": Check file size    OK")
        return True

def file_format(inFile):
    if inFile.endswith(".vcf.gz"):
        try:
            tmpFile = VariantFile(fileName, 'r')
            tmpFile.close()           
        except:
            logging.error(inFile + ": Not implemented in files compressed by bgzip.")
            return  False
        else:
            logging.info(inFile + ": Check vcf.gz format    OK")
            return  True
    else:
        logging.error(inFile + ": Not in .vcf.gz format.")
        return False

def chr_format(vcfFile, chrList):
    chr_list = [rec.chrom for rec in vcfFile.fetch()]
    chr_unique = list(set(chr_list))
    
    if len(chr_unique) == 1:
        logging.info(vcfFile.filename.decode() + ": Check chr spreated    OK")
    else:
        logging.error(vcfFile.filename.decode() + ": Include more than one chromsome.")
        return False

    if chr_unique[0][:3] == "chr":
        logging.info(vcfFile.filename.decode() + ": Check chr prefix    OK")
    else:
        logging.error(vcfFile.filename.decode() + ": #CHROM cloumn does not meet the required format(chr + number).")
        return False
    if int(chr_unique[0][3:]) in [i for i in range(23)]:
       pass
    else:
        logging.warn(vcfFile.filename.decode()+ ": ChinaMAP imputation server only support chr1-22. This file should be excluded.")
        return False
        
    if chr_unique[0] in chrList:
        logging.warn(vcfFile.filename.decode() + ": The chromosome already exists in your file list.")
    else:
        chrList.append(chr_unique[0])

    return chrList

def position_sort(vcfFile):
    pos_list = [rec.pos for rec in vcfFile.fetch()]
    pos_list_sort = pos_list.sort()
    if pos_list == pos_list_sort:
        logging.info(vcfFile.filename.decode() + ": Check pos sorted    OK")
    else:
        logging.error(vcfFile.filename.decode() + ": The file is not sorted.")
        return False

def ref_hg38(chrLength, vcfFile, chrList):
    try:
        length = vcfFile.header.contigs.get(chrList[-1]).length
    except:
        logging.warn(vcfFile.filename.decode() + ": Please make sure this file is from Hg38.")
        return False
    else:
        if length == chrLength[chrList[-1]]:
            logging.info(vcfFile.filename.decode() + ": Check hg38    OK")
            return True
        else:
            logging.error(vcfFile.filename.decode() + ": Chromsome length doesn't match Hg38 reference.")
            return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check input files for ChinaMAP imputation server.')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', type=str, default = None, required=True, help="Input .vcf.gz files path list.")
    args = parser.parse_args()
    inFileList = args.input
    
    LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename='ChinaMAP_checkVCF.log', level=logging.DEBUG, format=LOG_FORMAT)

    chrLength={"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718 ,"chr15": 101991189, "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167, "chr21": 46709983, "chr22": 50818468}

    fileList = []
    with open(inFileList, "r") as f:
        for line in f.readlines():
            fileList.append(line.strip('\n'))
    
    for fileName in fileList:
        if get_FileSize(fileName):
            pass
        else:
            fileList.remove(fileName)

        if file_format(fileName):
            pass
        else:
            fileList.remove(fileName)
    
    chrList = []

    for fileName in fileList:
        vcfFile = VariantFile(fileName, 'r')
        tmpChrList = chr_format(vcfFile, chrList)
        if tmpChrList:
            chrList = tmpChrList
        else:
            continue

        if ref_hg38(chrLength, vcfFile, chrList):
            pass
    
        vcfFile.close()