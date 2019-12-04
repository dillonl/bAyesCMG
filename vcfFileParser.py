import sys
import gzip
import pickle
import variant
from cyvcf2 import VCF
from enum import Enum

class VUSVCF:

    vepAnnotationKey = {'Location': 0, 'Allele': 1, 'SYMBOL': 2, 'IMPACT': 3, 'Consequence': 4, 'Protein_position': 5, 'Amino_acids': 6, 'Existing_variation': 7, 'IND': 8, 'ZYG': 9, 'ExACpLI': 10, 'REVEL': 11, 'DOMAINS': 12, 'CSN': 13, 'PUBMED': 14}

    # def __init__(self, families, vcfFilePath, priorProbability, oddsPathogenicity, exponent, gnomadDB, bs1AF, bs2AF):
    def __init__(self, families, vcfFilePath, priorProbability, oddsPathogenicity, exponent, gnomAD_AF_Threshold, REVEL_Threshold, outputVCF):
        scoresMap = { 'minPathogenicScore': 9 , 'minLikelyPathogenicScore' : 5 , 'maxLikelyBenignScore': -4 , 'maxBenignScore': -8 }
        self.outputVCF = outputVCF
        self.vcfFilePath = vcfFilePath
        self.families = families
        self.sampleIdxs = {}
        self.priorProbability = priorProbability
        self.oddsPathogenicity = oddsPathogenicity
        self.exponent = exponent
        self.gnomAD_AF_Threshold = gnomAD_AF_Threshold
        self.REVEL_Threshold = REVEL_Threshold


    def getSampleInfo(self, sampleName, formatField, vcfLineSplit):
        sampleInfoSplit = vcfLineSplit[self.sampleIdxs[sampleName]].split(':')
        sampleInfo = {}
        for i in range(0, len(formatField)):
            sampleInfo[formatFieldSplit[i]] = sampleInfoSplit[i]
        if '1' not in sampleInfo["GT"]:
            return None
        return

    def getCSQ(self, key, csq):
        if key not in self.vepAnnotationKey:
            return None
        else:
            return csq[self.vepAnnotationKey[key]]

    def processVariants(self):
        samples = {}
        for familyID in self.families:
            for sampleID in self.families[familyID].samples:
                samples[sampleID] = self.families[familyID].samples[sampleID]


        cyVCF = VCF(self.vcfFilePath)
        cyVCF.add_info_to_header({"ID": "Evidence_Codes", "Number": "1", "Type": "String", "Description": "All ACMG evidence codes that apply to this variant"})
        cyVCF.add_info_to_header({"ID": "Posterior_Pathogenic_Probability", "Number": "1", "Type": "String", "Description": "Posterior Pathogenic Probability"})
        for v in cyVCF:
            var = variant.Variant(v, self.gnomAD_AF_Threshold, self.REVEL_Threshold)
            posterior = self.getPosterior(var)
            v.INFO["Evidence_Codes"] = var.getEvidenceCodesString()
            v.INFO["Posterior_Pathogenic_Probability"] = str(self.getPosterior(var))
            # self.outputVCF.write(str(v))
            print(str(v))
            exit(0)

    def getPosterior(self, variant):
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
