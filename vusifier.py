import pandas as pd
import argparse
import vcfFileParser
import pedFileParser

def main():
    parser = argparse.ArgumentParser(description="Naive Bayes framework for ACMG/AMP variant classification")
    parser.add_argument('-v', action="store", dest="vcfFilePath", help="path to the VCF file", required=True)
    parser.add_argument('-p', action="store", dest="pedFilePath", help="path to the PED file", required=True)
    results = parser.parse_args()

    acmg_evidences = pd.read_csv('data/ACMG_evidences.csv')

    # variants = vcfFileParser.extractVariants(results.vcfFilePath)
    families = pedFileParser.parserPedFile(results.pedFilePath)
    vcf = vcfFileParser.VCF(families, results.vcfFilePath, acmg_evidences)
    vcf.processVariants()
    # print(variants)
    # print(families)

if __name__ == "__main__":
    main()
