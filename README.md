# Imputation-server

## ChinaMAP_checkVCF.py
ChinaMAP_checkVCF.py is a tool written in Python to check input VCF files before submit ChinaMAP imputation server task. It checks if your input files meets the requirements of the ChinaMAP imputation server.

### Requirments

* pysam
* python3.X

### Example

#### Recommend

```bash
python3 ChinaMAP_checkVCF.py -i -i chr2.vcf.gz chr3.vcf.gz -r mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz
```
#### If you don't want to check the sites percentage that are not in ChinaMAP reference panel

```bash
python3 ChinaMAP_checkVCF.py -i chr2.vcf.gz chr3.vcf.gz
```

### Result
* ChinaMAP_checkVCF.log

### Notice

The mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz file must has index file.
If the check fails, please go to http://www.mbiobank.com/imputation/help/ to find a solution.

## ChinaMAP_filterVCF.py
ChinaMAP_filterVCF.py is a tool written in Python to filter input VCF files before submit ChinaMAP imputation server task. It remove monomorphic sites and sites that are not in ChinaMAP reference panel from your input files.

### Requirments

* pysam
* python3.X
* mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz (You can download this file from http://www.mbiobank.com/download/)

### Example

#### If you want to remove sites that are not in ChinMAP reference panel

```bash
python3 ChinaMAP_filterVCF.py -i chr2.vcf.gz chr3.vcf.gz -d ./output/ -r mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz
```

#### If you want to remove monomorphic sites

```bash
python3 ChinaMAP_filterVCF.py -i chr2.vcf.gz chr3.vcf.gz -d ./output/ -e
```

#### If you want to remove both sites that not in ChinMAP reference panel and monomorphic sites

```bash
python3 ChinaMAP_filterVCF.py -i chr2.vcf.gz chr3.vcf.gz -d ./output/ -r mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz -e
```

#### Notice

The mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz file must has index file.
Please use at least one of -r or -e

### Result

* *_filtered.vcf.gz: Filtered ChinaMAP imputation server input files.
* *_exclude.vcf.gz: All sites that are not in ChinaMAP reference panel.
* *_monomorphic.vcf.gz: All monomorphic sites.