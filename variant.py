from cyvcf2 import VCF

class Variant:

    def __init__(self, variant, families, gnomAD_AF_Threshold, REVEL_Threshold, CSQList, matchingClinVarVariants):
        self.matchingClinVarVariants = matchingClinVarVariants
        self.variant = variant
        self.evidenceCodes = {"PVS1": 0, "PS1": 0, "PS2": 0, "PS3": 0, "PS4": 0, "PM1": 0, "PM2": 0, "PM3": 0, "PM4": 0, "PM5": 0, "PM6": 0, "PP1": 0, "PP2": 0, "PP3": 0, "PP4": 0, "PP5": 0, "BA1": 0, "BS1": 0, "BS2": 0, "BS3": 0, "BS4": 0, "BP1": 0, "BP2": 0, "BP3": 0, "BP4": 0, "BP5": 0, "BP6": 0, "BP7": 0}
        self.CSQList = CSQList
        self.CSQDict = self.getCSQDict(self.variant)
        '''
        self.CSQDict = {v:[] for i, v in enumerate(CSQList)} # produces a dictionary of empty arrays
        if self.variant.INFO.get("CSQ") != None:
            rawCSQ = [x.split('|') for x in self.variant.INFO.get("CSQ").split(',')]
            for csqIdx, csqKey in enumerate(CSQList):
                for csqList in rawCSQ:
                    self.CSQDict[csqKey].append(csqList[csqIdx])
        '''

        self.position = self.variant.POS
        self.gnomad_popmax_af = 0
        self.PStCounts = 0
        self.PMCounts = 0
        self.PSuCounts = 0
        self.PVstCounts = 0
        self.BStCounts = 0
        self.BSuCounts = 0
        self.BAsCounts = 0
        self.gnomAD_AF_Threshold = gnomAD_AF_Threshold
        self.REVEL_Threshold = REVEL_Threshold
        self.families = families
        self.parseLine()

    def parseLine(self):
        self.populateCSQ()
        self.processPVstCounts()
        self.processPStCounts()
        self.processPMCounts()
        self.processPSuCounts()

        self.processBAsCounts()
        self.processBStCounts()
        self.processBSuCounts()

    def getCSQDict(self, variant):
        CSQDict = {v:[] for i, v in enumerate(self.CSQList)} # produces a dictionary of empty arrays
        if variant.INFO.get("CSQ") != None:
            rawCSQ = [x.split('|') for x in variant.INFO.get("CSQ").split(',')]
            for csqIdx, csqKey in enumerate(self.CSQList):
                for csqList in rawCSQ:
                    CSQDict[csqKey].append(csqList[csqIdx])
        return CSQDict

    def processPVstCounts(self):
        # very strong PVS1
        highCount = self.CSQDict['IMPACT'].count("HIGH")
        if highCount > 0:
            motherGeno = self.variant.genotypes[self.families.getMotherIDFromVCFIdx()]
            fatherGeno = self.variant.genotypes[self.families.getFatherIDFromVCFIdx()]
            individualGeno = self.variant.genotypes[self.families.getIndividualIDFromVCFIdx()]
            if ((motherGeno[0] == 0 and motherGeno[1] == 0 and
                 fatherGeno[0] == 0 and fatherGeno[1] == 0) and
                ((individualGeno[0] == 0 and individualGeno[1] == 1) or (individualGeno[0] == 1 and individualGeno[1] == 0))):
                self.PStCounts = 1
                self.evidenceCodes["PVS1"] = 1
        else:
            self.evidenceCodes["PVS1"] = -1

    def processPStCounts(self):
        # strong PS1-PS4
        # Checking PS1
        self.evidenceCodes["PS1"] = 0
        if len(self.matchingClinVarVariants) > 0:
            for clinVarVariant in self.matchingClinVarVariants:
                if clinVarVariant.INFO.get('CLNSIG') is not None and clinVarVariant.INFO.get('CLNSIG').lower() == "pathogenic":
                    clinvarCSQDict = self.getCSQDict(clinVarVariant)
                    if clinvarCSQDict['Protein_position'] == self.CSQDict['Protein_position'] and clinvarCSQDict['Amino_acids'] == self.CSQDict['Amino_acids']:
                        self.evidenceCodes["PS1"] = 1
            if self.evidenceCodes["PS1"] != 1:
                self.evidenceCodes["PS1"] = -1
        # Checking PS2
        motherGeno = self.variant.genotypes[self.families.getMotherIDFromVCFIdx()]
        fatherGeno = self.variant.genotypes[self.families.getFatherIDFromVCFIdx()]
        individualGeno = self.variant.genotypes[self.families.getIndividualIDFromVCFIdx()]
        if ((motherGeno[0] == 0 and motherGeno[1] == 0 and
             fatherGeno[0] == 0 and fatherGeno[1] == 0) and
            ((individualGeno[0] == 0 and individualGeno[1] == 1) or (individualGeno[0] == 1 and individualGeno[1] == 0))):
            self.evidenceCodes["PS2"] = 1
        else:
            self.evidenceCodes["PS2"] = 0
        # Checking PS3
        self.evidenceCodes["PS3"] = 0
        # Checking PS4
        if (((individualGeno[0] == 0 and individualGeno[1] == 1) and (motherGeno[0] == 0 and motherGeno[1] == 0 and fatherGeno[0] == 0 and fatherGeno[1] == 0)) or ((individualGeno[0] == 1 and individualGeno[1] == 1) and (motherGeno[0] == 0 and motherGeno[1] == 1 and fatherGeno[0] == 0 and fatherGeno[1] == 1))):
            self.evidenceCodes["PS4"] = 1
        else:
            self.evidenceCodes["PS4"] = 0

    def processPMCounts(self):
        # moderate PM1-PM6
        # Checking PM1
        pfamCount = self.CSQDict['DOMAINS'].count("Pfam_domain")
        if pfamCount > 0:
            self.evidenceCodes["PM1"] = pfamCount
        else:
            self.evidenceCodes["PM1"] = -1
        # Checking PM2
        if self.gnomad_popmax_af < self.gnomAD_AF_Threshold:
            self.evidenceCodes["PM2"] = 1
        else:
            self.evidenceCodes["PM2"] = -1
        # Checking PM3
        self.evidenceCodes["PM3"] = 0
        # PM4 NA
        consequenceList = ['inframe_deletion', 'inframe_insertion', 'stop_lost']
        for consequence in self.CSQDict['Consequence']:
            if consequence in consequenceList:
                self.evidenceCodes["PM4"] = 1
        if self.evidenceCodes["PM4"] == 0:
            self.evidenceCodes["PM4"] = -1
        # PM5 NA
        self.evidenceCodes["PM5"] = 0
        if len(self.matchingClinVarVariants) > 0:
            for clinVarVariant in self.matchingClinVarVariants:
                if clinVarVariant.INFO.get('CLNSIG') is not None and clinVarVariant.INFO.get('CLNSIG').lower() == "pathogenic":
                    clinvarCSQDict = self.getCSQDict(clinVarVariant)
                    if clinvarCSQDict['Protein_position'] == self.CSQDict['Protein_position'] and clinvarCSQDict['Amino_acids'] != self.CSQDict['Amino_acids']:
                        print(self.variant)
                        self.evidenceCodes["PM5"] = 1
            if self.evidenceCodes["PM5"] != 1:
                self.evidenceCodes["PM5"] = 0
        # exit(0)

        # PM6 NA
        self.evidenceCodes["PM6"] = 0

    def processPSuCounts(self):
        # supporting PP1-PP5
        # Checking PP1
        self.evidenceCodes["PP1"] = 0
        # Checking PP2
        self.evidenceCodes["PP2"] = 0
        # Checking PP3
        self.evidenceCodes["PP3"] = 0
        for revel in self.CSQDict['REVEL']:
            if len(revel) > 0 and float(revel) > self.REVEL_Threshold:
                self.evidenceCodes["PP3"] = 1
        if self.evidenceCodes["PP3"] == 0:
            self.evidenceCodes["PP3"] = -1
        # Checking PP4
        self.evidenceCodes["PP4"] = 0
        # Checking PP5
        self.evidenceCodes["PP5"] = 0

    def processBAsCounts(self):
        # benign standalone BA1
        # Checking BA1
        if self.variant.INFO['gnomad_popmax_af'] > 0.05:
            self.evidenceCodes["BA1"] = 1
        else:
            self.evidenceCodes["BA1"] = -1

    def processBStCounts(self):
        motherGeno = self.variant.genotypes[self.families.getMotherIDFromVCFIdx()][:2]
        fatherGeno = self.variant.genotypes[self.families.getFatherIDFromVCFIdx()][:2]
        individualGeno = self.variant.genotypes[self.families.getIndividualIDFromVCFIdx()][:2]
        # strong BS1-BS4
        # Checking BS1
        if self.gnomad_popmax_af > self.gnomAD_AF_Threshold:
            self.evidenceCodes["BS1"] = 1
        else:
            self.evidenceCodes["BS1"] = -1
        # Checking BS2
        if ((individualGeno == motherGeno or individualGeno == fatherGeno) and (individualGeno == [0,1] or individualGeno == [1,1])):
            self.evidenceCodes["BS2"] = 1
        else:
            self.evidenceCodes["BS2"] = -1

        # Checking BS3
        self.evidenceCodes["BS3"] = 0
        # Checking BS4
        self.evidenceCodes["BS4"] = self.evidenceCodes["BS2"]

    def processBSuCounts(self):
        # supporting BP1-BP7
        # Checking BP1
        self.evidenceCodes["BP1"] = 0
        # Checking BP2
        self.evidenceCodes["BP2"] = 0
        # Checking BP3
        self.evidenceCodes["BP3"] = 0
        if (len(''.join(self.CSQDict['DOMAINS'])) == 0):
            consequenceList = ['inframe_deletion', 'inframe_insertion', 'stop_lost']
            for consequence in self.CSQDict['Consequence']:
                if consequence in consequenceList:
                    self.evidenceCodes["BP3"] = 1
        if self.evidenceCodes["BP3"] == 0:
            self.evidenceCodes["BP3"] = -1
        # Checking BP4
        self.evidenceCodes["BP4"] = -1
        for revel in self.CSQDict['REVEL']:
            if len(revel) > 0 and float(revel) < self.REVEL_Threshold:
                self.evidenceCodes["BP4"] = 1
                # Do we want to break here or keep counting?
        # Checking BP5
        self.evidenceCodes["BP5"] = 0
        # Checking BP6
        self.evidenceCodes["BP6"] = 0
        # Checking BP7
        self.evidenceCodes["BP7"] = self.CSQDict['Consequence'].count('synonymous_variant')# and self.CSQDict['Consequence'].count('splice_region')
        if self.evidenceCodes["BP7"] == 0 or (len(''.join(self.CSQDict['DOMAINS'])) > 0):
            self.evidenceCodes["BP7"] = -1

    def getEvidenceCodesString(self):
        codesOrder = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
        codes = []
        for code in codesOrder:
            codes.append(code + "=" + str(self.evidenceCodes[code]))
        return ';'.join(codes)

    def populateCSQ(self):
        #Location|Allele|SYMBOL|IMPACT|Consequence|Protein_position|Amino_acids|Existing_variation|IND|ZYG|ExACpLI|REVEL|DOMAINS|CSN|PUBMED
        idxs = ['Location', 'Allele', 'SYMBOL', 'IMPACT', 'Consequence', 'Protein_position', 'Amino_acids', 'Existing_variation', 'IND', 'ZYG', 'ExACpLI', 'REVEL', 'DOMAINS', 'CSN', 'PUBMED']
        csq = self.variant.INFO.get('CSQ')
        if csq is None: return
        csqSplit = csq.split("|")
        for i in range(len(csqSplit)):
            idx = i%len(idxs)
            key = idxs[idx]
        # print(self.variant.INFO.get('gnomad_popmax_af'))
        # print(self.variant)
        self.gnomad_popmax_af = float(self.variant.INFO.get('gnomad_popmax_af'))

    def getPathogenicSupportingCounts(self):
        return max(0, self.evidenceCodes["PP1"]) + max(0, self.evidenceCodes["PP2"]) + max(0, self.evidenceCodes["PP3"]) + max(0, self.evidenceCodes["PP4"]) + max(0, self.evidenceCodes["PP5"])

    def getPathogenicModerateCounts(self):
        return max(0, self.evidenceCodes["PM1"]) + max(0, self.evidenceCodes["PM2"]) + max(0, self.evidenceCodes["PM3"]) + max(0, self.evidenceCodes["PM4"]) + max(0, self.evidenceCodes["PM5"]) + max(0, self.evidenceCodes["PM6"])

    def getPathogenicStrongCounts(self):
        return max(0, self.evidenceCodes["PS1"]) + max(0, self.evidenceCodes["PS2"]) + max(0, self.evidenceCodes["PS3"]) + max(0, self.evidenceCodes["PS4"])

    def getPathogenicVeryStrongCounts(self):
        return max(0, self.evidenceCodes["PVS1"])

    def getBenignSupportingCounts(self):
        return max(0, self.evidenceCodes["BP1"]) + max(0, self.evidenceCodes["BP2"]) + max(0, self.evidenceCodes["BP3"]) + max(0, self.evidenceCodes["BP4"]) + max(0, self.evidenceCodes["BP5"]) + max(0, self.evidenceCodes["BP6"]) + max(0, self.evidenceCodes["BP7"])

    def getBenignModerateCounts(self):
        return 0

    def getBenignStrongCounts(self):
        return max(0, self.evidenceCodes["BS1"]) + max(0, self.evidenceCodes["BS2"]) + max(0, self.evidenceCodes["BS3"]) + max(0, self.evidenceCodes["BS4"])

    def getBenignVeryStrongCounts(self):
        return max(0, self.evidenceCodes["BA1"])
