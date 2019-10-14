import gzip

class VCF:

    def __init__(self, families, vcfFilePath):
        self.vcfFilePath = vcfFilePath
        self.families = families
        self.sampleIdxs = {}

    def getSampleInfo(self, sampleName, formatField, vcfLineSplit):
        sampleInfoSplit = vcfLineSplit[self.sampleIdxs[sampleName]].split(':')
        sampleInfo = {}
        for i in range(0, len(formatField)):
            sampleInfo[formatFieldSplit[i]] = sampleInfoSplit[i]
        if '1' not in sampleInfo["GT"]:
            return None
        return

    def processVariants(self):
        samples = {}
        for familyID in self.families:
            for sampleID in self.families[familyID].samples:
                samples[sampleID] = self.families[familyID].samples[sampleID]

        for line in gzip.open(self.vcfFilePath, 'rt'):
            line = line.replace('\n', '')
            lineSplit = line.split('\t')
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    counter = 9
                    for sampleName in lineSplit[9:]:
                        self.sampleIdxs[counter] = sampleName
                        counter += 1
                continue
            sampleLabelsSplit = lineSplit[-1*len(samples) - 1].split(":")
            for sampleIdx in range(9, 9+len(samples)):
                sampleSplit = lineSplit[sampleIdx].split(":")
                sampleDict = {sampleLabelsSplit[i] : sampleSplit[i] for i in range(len(sampleSplit)) }
                print(sampleDict)
                # AT THIS POINT YOU"LL WANT TO GET THE EVIDENCE CODES FOR EACH AND PASS THEM TO
                # A COOL FUNCTION THAT CALCULATES THE POSTERIOR PROBABILITY AND THE ODDS_PATH
            exit(0)
            '''
            info = {f.split("=")[0]: f.split("=")[1] if "=" in f else f for f in lineSplit[7].split(';')}
            if 'CSQ' not in info:
                continue
            csqValues = info["CSQ"].split(",")
            print(len(self.vep_keys), len(csqValues))
            print(info["CSQ"])
            csq = {self.vep_keys[i]: csqValues[i] for i in range(len(csqValues))}
            print(csq)
            exit(0)


            formatField = lineSplit[8].split(':')
            variantSamples = {}
            # for sample in samples:
                # variantSamples[sample] = {key:val for }
            '''
