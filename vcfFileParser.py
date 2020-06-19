import sys
import gzip
import pickle
import variant
from cyvcf2 import VCF
from enum import Enum

from pprint import pprint

class VUSVCF:

    vepAnnotationKey = {'Location': 0, 'Allele': 1, 'SYMBOL': 2, 'IMPACT': 3, 'Consequence': 4, 'Protein_position': 5, 'Amino_acids': 6, 'Existing_variation': 7, 'IND': 8, 'ZYG': 9, 'ExACpLI': 10, 'REVEL': 11, 'DOMAINS': 12, 'CSN': 13, 'PUBMED': 14}

    def __init__(self, families, vcfFilePath, priorProbability, oddsPathogenicity, exponent, gnomAD_AF_Threshold, REVEL_Threshold, outputVCF, clinVarData):
        scoresMap = { 'minPathogenicScore': 9 , 'minLikelyPathogenicScore' : 5 , 'maxLikelyBenignScore': -4 , 'maxBenignScore': -8 }
        self.clinVarData = clinVarData
        self.outputVCF = outputVCF
        self.vcfFilePath = vcfFilePath
        self.families = families
        self.sampleIdxs = {}
        self.priorProbability = priorProbability
        self.oddsPathogenicity = oddsPathogenicity
        self.exponent = exponent
        self.gnomAD_AF_Threshold = gnomAD_AF_Threshold
        self.REVEL_Threshold = REVEL_Threshold

        self.evidenceExponentSupporting = self.exponent ** (-3) # found in the spreadsheet
        self.evidenceExponentModerate = self.exponent ** (-2) # found in the spreadsheet
        self.evidenceExponentStrong = self.exponent ** (-1) # found in the spreadsheet
        self.evidenceExponentVeryStrong = self.exponent ** (0) # found in the spreadsheet

        self.oddsPathSupporting = oddsPathogenicity ** self.evidenceExponentSupporting
        self.oddsPathModerate = oddsPathogenicity ** self.evidenceExponentModerate
        self.oddsPathStrong = oddsPathogenicity ** self.evidenceExponentStrong
        self.oddsPathVeryStrong = oddsPathogenicity ** self.evidenceExponentVeryStrong


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
        cyVCF = VCF(self.vcfFilePath)
        self.families.setSampleIdxs(cyVCF.samples)
        getCSQList = self.getCSQList(cyVCF.raw_header)
        cyVCF.add_info_to_header({"ID": "Evidence_Codes", "Number": "1", "Type": "String", "Description": "All ACMG evidence codes that apply to this variant"})
        cyVCF.add_info_to_header({"ID": "Posterior_Pathogenic_Probability", "Number": "1", "Type": "String", "Description": "Posterior Pathogenic Probability"})
        self.outputVCF.write(cyVCF.raw_header)
        for v in cyVCF:
            matchingClinVarVariants = []
            for alt in v.ALT:
                key = "%s:%s:%s:%s" % (v.CHROM, v.POS, v.REF, alt)
                if key in self.clinVarData:
                    matchingClinVarVariants.append(self.clinVarData[key])
            var = variant.Variant(v, self.families, self.gnomAD_AF_Threshold, self.REVEL_Threshold, getCSQList, matchingClinVarVariants)
            if not var.printVariant:
                continue
            posterior = self.getPosterior(var)
            v.INFO["Evidence_Codes"] = var.getEvidenceCodesString()
            v.INFO["Posterior_Pathogenic_Probability"] = str(format(self.getPosterior(var), '.3f'))
            self.outputVCF.write(str(v))

    def getCSQList(self, rawHeader):
        headerSplit1 = rawHeader.split('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: ')
        headerSplit2 = headerSplit1[1].split('">')
        return headerSplit2[0].split('|')

    def getPosterior(self, variant):
        pathogenicOddsPathSupporting = self.oddsPathSupporting ** variant.getPathogenicSupportingCounts()
        pathogenicOddsPathModerate = self.oddsPathModerate ** variant.getPathogenicModerateCounts()
        pathogenicOddsPathStrong = self.oddsPathStrong ** variant.getPathogenicStrongCounts()
        pathogenicOddsPathVeryStrong = self.oddsPathVeryStrong ** variant.getPathogenicVeryStrongCounts()

        benignOddsPathSupporting = self.oddsPathSupporting ** ((-1)*variant.getBenignSupportingCounts())
        benignOddsPathModerate = self.oddsPathModerate ** ((-1)*variant.getBenignModerateCounts())
        benignOddsPathStrong = self.oddsPathStrong ** ((-1)*variant.getBenignStrongCounts())
        benignOddsPathVeryStrong = self.oddsPathVeryStrong ** ((-1)*variant.getBenignVeryStrongCounts())

        combinedOdds = pathogenicOddsPathSupporting * pathogenicOddsPathModerate * pathogenicOddsPathStrong * pathogenicOddsPathVeryStrong * benignOddsPathSupporting * benignOddsPathModerate * benignOddsPathStrong * benignOddsPathVeryStrong

        posteriorNum = (float(combinedOdds) * float(self.priorProbability))
        posteriorDenom = (((float(combinedOdds)-1)*float(self.priorProbability))+float(1))
        posterior = posteriorNum/posteriorDenom

        # print(benignOddsPathStrong, variant.getBenignStrongCounts())
        # print(benignOddsPathVeryStrong, variant.getBenignVeryStrongCounts())

        # print(pathogenicOddsPathSupporting, pathogenicOddsPathModerate, pathogenicOddsPathStrong, pathogenicOddsPathVeryStrong, benignOddsPathSupporting, benignOddsPathModerate, benignOddsPathStrong, benignOddsPathVeryStrong)
        # print(combinedOdds, posterior)

        return posterior
