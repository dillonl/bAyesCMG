from cyvcf2 import VCF
import sys

if __name__ == '__main__':
    vcfFile = sys.argv[1]
    cyVCF = VCF(vcfFile)

    if 'ID=CSQ' in cyVCF.raw_header and 'ID=CLNSIG' in cyVCF.raw_header:
        exit(0)
    else:
        exit(1)
