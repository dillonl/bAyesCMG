import sys
import gzip
# sys.path.append('externals/CharGer')
sys.path.append('externals/CharGer/charger')
import jimmer
from charger import charger
import imp
imp.reload(charger)
from enum import Enum

class CODES(Enum):
    #benign
    BP1 = 1
    BP2 = 2
    BP3 = 3
    BP4 = 4
    BP5 = 5
    BP6 = 6
    BP7 = 7
    BS1 = 8
    BS2 = 9
    BS3 = 10
    BS4 = 11
    BA1 = 12
    #pathogenic
    PP1 = 13
    PP2 = 14
    PP3 = 15
    PM1 = 16
    PM2 = 17
    PM3 = 18
    PM4 = 19
    PM5 = 20
    PM6 = 21
    PS1 = 22
    PS2 = 23
    PS3 = 24
    PS4 = 25
    PVS1 = 26

class VCF:

    def __init__(self, families, vcfFilePath, priorProbability, oddsPathogenicity, exponent):
        self.vcfFilePath = vcfFilePath
        self.families = families
        self.sampleIdxs = {}
        self.priorProbability = priorProbability
        self.oddsPathogenicity = oddsPathogenicity
        self.exponent = exponent


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



    def getPosterior(variant):
        countPVst = variant.getPVstCounts()
        countPSt = variant.getPStCounts()
        countPM = variant.getPMCounts()
        countPSu = variant.getPSuCounts()
        oddsPathEquation = self.oddsPathogenicity ** ((countPSu / (self.exponent ** 3)) + (countPM / (self.exponent ** 2)) + (countPSt / (self.exponent ** 1)) + (countPVst / 1))

        countBSt = variant.getBStCounts()
        countBSu = variant.getBSuCounts()
        oddsBenignEquation = self.oddsPathogenicity ** ((countBSu / (self.exponent ** 3)) + (countBSt / (self.exponent ** 2)))

        combinedOddsPath = oddsPathEquation * oddsBenignEquation

        posterior = (combinedOddsPath * self.priorProbability) / ((combinedOddsPath - 1) * (self.priorProbability + 1))
        return posterior
