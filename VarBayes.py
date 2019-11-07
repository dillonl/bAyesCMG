import sys
import pandas as pd
import argparse
import vcfFileParser
import pedFileParser
import gzip
import os.path
from ftplib import FTP
from dateutil import parser
from datetime import datetime, timedelta
from pytz import timezone
sys.path.append('externals/CharGer')
from charger import charger
reload(charger)

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

def generateCharGerVCF():
    self.charger = jimmer.charger()
    # self.charger = charger.charger()
    [ vepDone , preVEP , exacDone , clinvarDone ] = self.charger.getInputData( vcf=vcfFilePath )
    self.charger.getExternalData( clinvar=True , exac=True , vep=True )
    self.charger.PVS1( )
    self.charger.PS1( )
    self.charger.PS2( )
    self.charger.PS3( )
    self.charger.PS4( )
    # self.charger.PM1( recurrenceThreshold , hotspot3d=clustersFile )
    # self.charger.PM2( rareAF )
    self.charger.PM3( )
    self.charger.PM4( )
    self.charger.PM5( )
    self.charger.PM6( )
    self.charger.PP1( )
    self.charger.PP2( )
    # self.charger.PP3( minimumEvidence )
    self.charger.PP4( )
    self.charger.PP5( )

    # self.charger.BA1( commonAF )
    self.charger.BS1( )
    self.charger.BS2( )
    self.charger.BS3( )
    self.charger.BS4( )
    self.charger.BP1( )
    self.charger.BP2( )
    self.charger.BP3( )
    # self.charger.BP4( minimumEvidence )
    self.charger.BP5( )
    self.charger.BP6( )
    self.charger.BP7( )

    self.charger.PSC1( )
    self.charger.PMC1( )
    self.charger.PPC1( )
    self.charger.PPC2( )

    self.charger.BSC1( )
    self.charger.BMC1( )

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
    vcf = vcfFileParser.VCF(families, results.vcfFilePath, results.priorProbability, results.oddsPathogenicity, results.exponent)
    vcf.processVariants()
    # print(variants)
    # print(families)

if __name__ == "__main__":
    main()
