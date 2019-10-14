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

def downloadClinVarIfNeeded(clinVarFileName):
    shouldDownloadClinVar = True
    if os.path.isfile(clinVarFileName):
        clinVarDate = None
        for l in gzip.open(clinVarFileName, 'rt'):
            if l.startswith('##fileDate='):
                l = l.replace('\n', '')
                clinVarDate = l.replace('##fileDate=', '')
                ftp = FTP('ftp.ncbi.nlm.nih.gov')
                ftp.login()
                ftp.cwd('/pub/clinvar/vcf_GRCh37/')
                timestamp = ftp.sendcmd('MDTM clinvar.vcf.gz')[4:].strip()
                uploadtime = parser.parse(timestamp)

                clinVarDateObj = datetime.strptime(clinVarDate, '%Y-%m-%d')
                diffTime = clinVarDateObj - uploadtime
                shouldDownloadClinVar = abs(diffTime.total_seconds()) > abs(timedelta(hours=48).total_seconds()) # since we don't have the server's timezone we don't know when it was uploaded so we need to add a fudgefactor 48 hours
                ftp.quit()
                break

    if shouldDownloadClinVar:
        clinVarURL = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
        # response = urllib.request.urlretrieve(clinVarURL, clinVarFileName)


def main():
    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    parser = argparse.ArgumentParser(description="Naive Bayes framework for ACMG/AMP variant classification")
    parser.add_argument('-v', action="store", dest="vcfFilePath", help="path to the VCF file", required=True)
    parser.add_argument('-p', action="store", dest="pedFilePath", help="path to the PED file", required=True)
    parser.add_argument('-a', action="store_true", dest="addAnnotations", help="use CharGer to generate vcf annotations", required=False)
    results = parser.parse_args()

    print(results.addAnnotations)
    exit(0)


    clinVarFileName = scriptPath + '/data/clinvar.vcf.gz'
    downloadClinVarIfNeeded(clinVarFileName)

    acmg_evidences = pd.read_csv(scriptPath + '/data/ACMG_evidences.csv')

    # variants = vcfFileParser.extractVariants(results.vcfFilePath)
    families = pedFileParser.parserPedFile(results.pedFilePath)
    vcf = vcfFileParser.VCF(families, results.vcfFilePath, acmg_evidences)
    vcf.processVariants()
    # print(variants)
    # print(families)

if __name__ == "__main__":
    main()
