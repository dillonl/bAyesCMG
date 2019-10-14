import pandas as pd
import argparse
import vcfFileParser
import pedFileParser
# import urllib.request
import gzip
import os.path
import sys
from ftplib import FTP
from dateutil import parser
from datetime import datetime, timedelta
from pytz import timezone
sys.path.append('externals/CharGer')
from charger import charger

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

def main():
    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    parser = argparse.ArgumentParser(description="Naive Bayes framework for ACMG/AMP variant classification")
    parser.add_argument('-v', action="store", dest="vcfFilePath", help="path to the VCF file", required=True)
    parser.add_argument('-f', action="store", dest="pedFilePath", help="path to the PED file", required=True)
    parser.add_argument('-p', action="store", dest="priorProbability", help='The prior probability (default value = 0.1)', required=False, default=0.1)
    parser.add_argument('-o', action="store", dest="oddsPathogenicity", help='The odds of pathogenicity for "Very Strong" (default value = 350)', required=False, default=350)
    parser.add_argument('-e', action="store", dest="exponent", help='The exponent that sets the strength of Supporting/ Moderate/ Strong/ compared to "Very Strong" (default value = 2.0)', required=False, default=2.0)
    parser.add_argument('-a', action="store_true", dest="addAnnotations", help="use CharGer to generate vcf annotations", required=False)
    results = parser.parse_args()
    if not validateCommandLineArgs(results):
        exit(0)

    families = pedFileParser.parserPedFile(results.pedFilePath)
    vcf = vcfFileParser.VCF(families, results.vcfFilePath)
    vcf.processVariants()
    # print(variants)
    # print(families)

if __name__ == "__main__":
    main()
