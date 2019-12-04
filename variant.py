from cyvcf2 import VCF

class Variant:

    def __init__(self, variant, gnomAD_AF_Threshold, REVEL_Threshold):
        self.variant = variant
        self.evidenceCodes = {}
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
        self.CSQ = {}
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

    def processPVstCounts(self):
        # very strong PVS1
        if "HIGH" in self.CSQ['Consequence']:
            self.PStCounts += 1
            self.evidenceCodes["PVS1"] = 1
        else:
            self.evidenceCodes["PVS1"] = -1

    def processPStCounts(self):
        # strong PS1-PS4
        # Checking PS1
        self.evidenceCodes["PS1"] = 0
        # Checking PS2
        self.evidenceCodes["PS2"] = 0
        # Checking PS3
        self.evidenceCodes["PS3"] = 0
        # Checking PS4
        self.evidenceCodes["PS4"] = 0

    def processPMCounts(self):
        # moderate PM1-PM6
        # Checking PM1
        self.evidenceCodes["PM1"] = -1
        for d in self.CSQ['DOMAINS']:
            if 'Pfam_domain' in d:
                self.evidenceCodes["PM1"] = 1
                print('PM1 Detected to learn more email: Matt Velinder <mvelinder@gmail.com>:', d)
        # Checking PM2
        if self.gnomad_popmax_af < self.gnomAD_AF_Threshold:
            self.evidenceCodes["PM2"] = 1
        else:
            self.evidenceCodes["PM2"] = -1
        # Checking PM3
        self.evidenceCodes["PM3"] = -1
        for revel in self.CSQ['REVEL']:
            if len(revel) > 0 and float(revel) > self.REVEL_Threshold:
                self.evidenceCodes["PM3"] = 1
                # Do we want to break here or keep counting?
        # PM4 NA
        self.evidenceCodes["PM4"] = 0
        # PM5 NA
        self.evidenceCodes["PM5"] = 0
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
        # Checking PP4
        self.evidenceCodes["PP4"] = 0
        # Checking PP5
        self.evidenceCodes["PP5"] = 0

    def processBAsCounts(self):
        # benign standalone BA1
        # Checking BA1
        if self.gnomad_popmax_af < 0.05:
            self.evidenceCodes["BA1"] = 1
        else:
            self.evidenceCodes["BA1"] = -1

    def processBStCounts(self):
        # strong BS1-BS4
        # Checking BS1
        if self.gnomad_popmax_af > self.gnomAD_AF_Threshold:
            self.evidenceCodes["BS1"] = 1
        else:
            self.evidenceCodes["BS1"] = -1
        # Checking BS2
        self.evidenceCodes["BS2"] = 0
        # Checking BS3
        self.evidenceCodes["BS3"] = 0
        # Checking BS4
        self.evidenceCodes["BS4"] = 0

    def processBSuCounts(self):
        # supporting BP1-BP7
        # Checking BP1
        self.evidenceCodes["BP1"] = 0
        # Checking BP2
        self.evidenceCodes["BP2"] = 0
        # Checking BP3
        self.evidenceCodes["BP3"] = 0
        # Checking BP4
        self.evidenceCodes["BP4"] = -1
        for revel in self.CSQ['REVEL']:
            if len(revel) > 0 and float(revel) < self.REVEL_Threshold:
                self.evidenceCodes["BP4"] = 1
                # Do we want to break here or keep counting?
        # Checking BP5
        self.evidenceCodes["BP5"] = 0
        # Checking BP6
        self.evidenceCodes["BP6"] = 0
        # Checking BP7
        self.evidenceCodes["BP7"] = 0

    def getEvidenceCodesString(self):
        codesOrder = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
        codes = []
        for code in codesOrder:
            codes.append(code + "=" + str(self.evidenceCodes[code]))
        return ';'.join(codes)


    def getPVstCounts(self):
        # very strong PVS1
        return max(0, self.evidenceCodes["PVS1"])

    def getPStCounts(self):
        # strong PS1-PS4
        return max(0, self.evidenceCodes["PS1"]) + max(0, self.evidenceCodes["PS2"]) + max(0, self.evidenceCodes["PS3"]) + max(0, self.evidenceCodes["PS4"])

    def getPMCounts(self):
        # supporting PM1-PM6
        return max(0, self.evidenceCodes["PM1"]) + max(0, self.evidenceCodes["PM2"]) + max(0, self.evidenceCodes["PM3"]) + max(0, self.evidenceCodes["PM4"]) + max(0, self.evidenceCodes["PM5"]) + max(0, self.evidenceCodes["PM6"])

    def getPSuCounts(self):
        # supporting PP1-PP5
        return max(0, self.evidenceCodes["PP1"]) + max(0, self.evidenceCodes["PP2"]) + max(0, self.evidenceCodes["PP3"]) + max(0, self.evidenceCodes["PP4"]) + max(0, self.evidenceCodes["PP5"])

    def getBAsCounts(self):
        # benign standalone BA1
        return max(0, self.evidenceCodes["BA1"])

    def getBStCounts(self):
        # strong BS1-BS4
        return max(0, self.evidenceCodes["BS1"]) + max(0, self.evidenceCodes["BS2"]) + max(0, self.evidenceCodes["BS3"]) + max(0, self.evidenceCodes["BS4"])

    def getBSuCounts(self):
        # supporting BP1-BP7
        return self.BSuCounts

    def populateCSQ(self):
        #Location|Allele|SYMBOL|IMPACT|Consequence|Protein_position|Amino_acids|Existing_variation|IND|ZYG|ExACpLI|REVEL|DOMAINS|CSN|PUBMED
        idxs = ['Location', 'Allele', 'SYMBOL', 'IMPACT', 'Consequence', 'Protein_position', 'Amino_acids', 'Existing_variation', 'IND', 'ZYG', 'ExACpLI', 'REVEL', 'DOMAINS', 'CSN', 'PUBMED']
        self.CSQ = {'Location': [], 'Allele': [], 'SYMBOL': [], 'IMPACT': [], 'Consequence': [], 'Protein_position': [], 'Amino_acids': [], 'Existing_variation': [], 'IND': [], 'ZYG': [], 'ExACpLI': [], 'REVEL': [], 'DOMAINS': [], 'CSN': [], 'PUBMED': []}
        csq = self.variant.INFO.get('CSQ')
        if csq is None: return
        csqSplit = csq.split("|")
        for i in range(len(csqSplit)):
            idx = i%len(idxs)
            key = idxs[idx]
            self.CSQ[key].append(csqSplit[idx])
        self.gnomad_popmax_af = float(self.variant.INFO.get('gnomad_popmax_af'))
