import gzip

class Family:

    def __init__(self, familyID):
        self.familyID = familyID
        self.paternalSample = None
        self.maternalSample = None
        self.childSample = None
        self.AllFamilySamplesDict = {}

    def setPaternalSample(self, paternalSample):
        self.AllFamilySamplesDict[paternalSample.individualID] = paternalSample
        self.paternalSample = paternalSample

    def setMaternalSample(self, maternalSample):
        self.AllFamilySamplesDict[maternalSample.individualID] = maternalSample
        self.maternalSample = maternalSample

    def getProband(self):
        affectedSamples = self.getAllAffectedIndividuals()
        for affSamp in affectedSamples:
            if affSamp.maternalSample != None and affSamp.paternalSample != None:
                return affSam
        return None

    def setChildSample(self, childSample):
        self.AllFamilySamplesDict[childSample.individualID] = childSample
        self.childSample = childSample

    def getAllAffectedIndividuals(self):
        affectedSamples = []
        for individualID, sample in self.AllFamilySamplesDict.items():
            if sample.affected:
                affectedSamples.append(sample)
        return affectedSamples

    def getAllUnaffectedIndividuals(self):
        unaffectedSamples = []
        for individualID, sample in self.AllFamilySamplesDict.items():
            if not sample.affected:
                unaffectedSamples.append(sample)
        return unaffectedSamples

    def getAllFamilySamples(self):
        familySamples = []
        for individualID, sample in self.AllFamilySamplesDict.items():
            familySamples.append(sample)
        return familySamples

    def setSampleIdxs(self, samplesFromVCF):
        for i, sampleID in enumerate(samplesFromVCF):
            self.AllFamilySamplesDict[sampleID].setSampleVCFIdx(i)

    def getMaternalSample(self):
        return self.maternalSample
    def getPaternalSample(self):
        return self.paternalSample
    def getChildSample(self):
        return self.childSample

class Sample:

    def __init__(self, individualID, paternalID, maternalID, sex, phenotype):
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
        self.affected = (self.phenotype == 2)

    def setSampleVCFIdx(self, idx):
        self.sampleIdx = idx

def parserPedFile(pedFilePath):
    families = {}
    for line in open(pedFilePath):
        if 'Kindred_ID' in line:
            continue
        line = line.replace('\n', '')
        lineSplit = line.split()
        if lineSplit[0] not in families:
            families[lineSplit[0]] = Family(lineSplit[0])
        sample = Sample(lineSplit[1], lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
        if sample.paternalID == '0' and sample.maternalID == '0':
            if sample.sex == 1:
                families[lineSplit[0]].setPaternalSample(sample)
            else:
                families[lineSplit[0]].setMaternalSample(sample)
        else:
            families[lineSplit[0]].setChildSample(sample)
    for key in families:
        return families[key] #just return the first family
