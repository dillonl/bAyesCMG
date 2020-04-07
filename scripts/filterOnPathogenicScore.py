import sys
import os
import gzip
import argparse

def main():
    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    parser = argparse.ArgumentParser(description="A tool to filter bAyesCMG VCFs on posterior probability.")
    parser.add_argument('-v', action="store", dest="vcfFilePath", help="path to the bAyesCMG VCF file", required=True)
    parser.add_argument('-t', action="store", dest="threshold", help='The posterior probability threshold, everything less than this will be filtered out (default value = 0.75)', required=True, default=0.75, type=float)
    results = parser.parse_args()
    vcfFilePath = results.vcfFilePath
    threshold = results.threshold
    if threshold < 0 or threshold > 1:
        print("Threshold needs to be between 0 and 1 (inclusive)")
        exit(0)
    if not os.path.isfile(vcfFilePath):
        print("You must provide a path to a valid VCF")
        exit(0)

    f = None
    if vcfFilePath.endswith('vcf'):
        f = open(vcfFilePath, 'r')
    elif vcfFilePath.endswith('gz'):
        f = gzip.open(vcfFilePath, 'rb')

    for line in f:
        line = line.replace('\n', '')
        if line.startswith('#'):
            print(line)
            continue
        lineSplit = line.split('\t')
        d = lineSplit[7].split('Posterior_Pathogenic_Probability=')
        if len(d) < 2:
            continue
        posterior = d[1]
        if ';' in posterior:
            posterior = posterior.split(';')[0]
        if float(posterior) >= threshold:
            print(line)

if __name__ == "__main__":
    main()
