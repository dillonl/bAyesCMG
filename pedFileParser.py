import gzip

class Family:

    def __init__(self, familyID, individualID, paternalID, maternalID):
        self.familyID = familyID
        self.individualID = individualID
        self.paternalID = paternalID
        self.maternalID = maternalID
        self.samples = {}
        self.IDtoIdx = {self.maternalID:0, self.paternalID:0, self.individualID:0}

    def getAllIndividuals(self):
        return self.samples

    def addIndividual(self, sample):
        if sample.individualID not in self.samples:
            self.samples[sample.individualID] = sample

    def getIndividual(self, individualID):
        return self.samples[individualID]

    def setSampleIdxs(self, samplesFromVCF):
        for i, sampleID in enumerate(samplesFromVCF):
            self.IDtoIdx[sampleID] = i

    def getMotherID(self):
        return self.maternalID
    def getFatherID(self):
        return self.paternalID
    def getIndividualID(self):
        return self.individualID

    def getMotherIDFromVCFIdx(self):
        return self.IDtoIdx[self.maternalID]
    def getFatherIDFromVCFIdx(self):
        return self.IDtoIdx[self.paternalID]
    def getIndividualIDFromVCFIdx(self):
        return self.IDtoIdx[self.individualID]

class Sample:

    def __init__(self, family, individualID, paternalID, maternalID, sex, phenotype):
        self.family = family
        self.individualID = individualID
        self.paternalID = paternalID
        self.maternalID = maternalID

        self.sampleIdx = 0
        try:
            self.sex = int(sex)
            self.phenotype = int(phenotype)
        except ValueError:
            print("PED file format is invalid. Program exiting...")
            exit(0)
        self.family.addIndividual(self)

def parserPedFile(pedFilePath):
    families = {}
    for line in open(pedFilePath):
        if line.startswith('#'):
            continue
        line = line.replace('\n', '')
        lineSplit = line.split()
        if lineSplit[0] not in families:
            families[lineSplit[0]] = Family(lineSplit[0], lineSplit[1], lineSplit[2], lineSplit[3])
        sample = Sample(families[lineSplit[0]], lineSplit[1], lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
    for key in families:
        return families[key] #just return the first family
'''
def parserPedFile(pedFilePath):
    families = {}
    for line in open(pedFilePath):
        if line.startswith('#'):
            continue
        line = line.replace('\n', '')
        lineSplit = line.split()
        if lineSplit[0] not in families:
            families[lineSplit[0]] = Family(lineSplit[0])
        sample = Sample(families[lineSplit[0]], lineSplit[1], lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
    return families
'''
