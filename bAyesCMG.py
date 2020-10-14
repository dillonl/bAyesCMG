import sys
import pandas as pd
import argparse
import vcfFileParser
import pedFileParser
import os.path
from cyvcf2 import VCF
from ftplib import FTP
from dateutil import parser
from datetime import datetime, timedelta
from pytz import timezone

def validateCommandLineArgs(results):
    if results.priorProbability != None:
        printPriorError = False
        try:
            priorProbability = float(results.priorProbability)
            printPriorError = not (0 <= priorProbability and priorProbability <= 1)
        except ValueError as err:
            printPriorError = True
        if printPriorError:
            print('Prior probability must be a float between 0 and 1 (inclusive)')
            return False

    oddsPathogenicity = 350
    if results.oddsPathogenicity != None:
        printOddsError = False
        try:
            priorProbability = int(results.oddsPathogenicity)
            printOddsError = not (0 <= oddsPathogenicity)
        except ValueError as err:
            printOddsError = True
        if printOddsError:
            print('The odds of pathogenicity must be a number greater than 0')
            return False

    if results.exponent != None:
        try:
            exponenet = float(results.exponent)
        except ValueError as err:
            print('exponent must be a numeric value')
            return False
    return True

def getClinVarData(clinVarPath):
    clinVarData = {}
    cyVCF = VCF(clinVarPath)
    for v in cyVCF:
        for alt in v.ALT:
            key = "%s:%s:%s:%s" % (v.CHROM, v.POS, v.REF, alt)
            clinVarData[key] = v
    return clinVarData


def main():

    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    parser = argparse.ArgumentParser(description="Naive Bayes framework for ACMG/AMP variant classification")
    parser.add_argument('-d', action="store", dest="outputVcfFilePath", help="path to the output VCF file", required=True)
    parser.add_argument('-v', action="store", dest="vcfFilePath", help="path to the VCF file", required=True)
    parser.add_argument('-f', action="store", dest="pedFilePath", help="path to the PED file", required=True)
    parser.add_argument('-c', action="store", dest="clinVar", help="Path to VEP annotated ClinVar file", required=True)
    parser.add_argument('-p', action="store", dest="priorProbability", help='The prior probability (default value = 0.1)', required=False, default=0.1, type=float)
    parser.add_argument('-o', action="store", dest="oddsPathogenicity", help='The odds of pathogenicity for "Very Strong" (default value = 350)', required=False, default=350, type=float)
    parser.add_argument('-e', action="store", dest="exponent", help='The exponent that sets the strength of Supporting/ Moderate/ Strong/ compared to "Very Strong" (default value = 2.0)', required=False, default=2.0, type=float)
    parser.add_argument('-a', action="store", dest="gnomAD_AF_Threshold", help='gnomAD_AF threshold [Ask Matt] (default value = 0.01)', required=False, default=0.01, type=float)
    parser.add_argument('-r', action="store", dest="REVEL_Threshold", help='REVEL threshold [Ask Matt] (default value = 0.6)', required=False, default=0.6, type=float)
    results = parser.parse_args()
    if not validateCommandLineArgs(results):
        exit(0)

    clinVarData = getClinVarData(results.clinVar)
    outputDirectory = results.outputVcfFilePath
    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)
    outputVCFFilePath = outputDirectory + "/bayescmg.vcf"
    outputVCFFile = open(outputVCFFilePath, 'w')
    families = pedFileParser.parserPedFile(results.pedFilePath)
    vcf = vcfFileParser.VUSVCF(families, results.vcfFilePath, results.priorProbability, results.oddsPathogenicity, results.exponent, results.gnomAD_AF_Threshold, results.REVEL_Threshold, outputVCFFile, clinVarData)
    vcf.processVariants()
    outputVCFFile.close()


if __name__ == "__main__":
    main()
