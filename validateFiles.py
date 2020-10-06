import sys, gzip

def getPEDSamples(pedFile):
    samples = []
    firstLine = True
    if pedFile.endswith('.gz'):
        f = gzip.open(pedFile, 'rb')
    else:
        f = open(pedFile, 'r')
    for l in f:
        if firstLine:
            firstLine = False
            continue
        samples.append(l.split()[1])
    return samples

def getVCFSamples(vcfFile):
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    samples = []
    for l in open(vcfFile, 'r'):
        if l.startswith('#CHROM'):
            for h in l.replace('\n', '').split():
                if h not in cols:
                    samples.append(h)
            return samples
        else:
            continue

def isPEDSamplesSubsetVCFSamples(pedSamples, vcfSamples):
    return all(x in vcfSamples for x in pedSamples)

def isVCFDecomposed(vcfFile):
    for l in open(vcfFile, 'r'):
        if l.startswith('#'):
            continue
        l = l.split('\t')
        if ',' in l[4]:
            return False
    return True

if __name__ == '__main__':
    vcfFile = sys.argv[1]
    pedFile = sys.argv[2]
    pedSamples = getPEDSamples(pedFile)
    vcfSamples = getVCFSamples(vcfFile)
    if not isPEDSamplesSubsetVCFSamples(pedSamples, vcfSamples):
        print('PED samples are not a subset of the VCF samples')
        print('PED samples', pedSamples)
        print('VCF samples', vcfSamples)
        exit(1)
    '''
    if not isVCFDecomposed(vcfSamples):
        print('VCF is not decomposed')
        exit(1)
    '''
    exit(0)
