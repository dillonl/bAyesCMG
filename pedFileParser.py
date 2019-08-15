import gzip

class Family:

    def __init__(self, familyID):
        self.familyID = familyID
        self.samples = {}

    def getAllIndividuals(self):
        return self.samples

    def addIndividual(self, sample):
        if sample.individualID not in self.samples:
            self.samples[sample.individualID] = sample

    def getIndividual(self, individualID):
        return self.samples[individualID]

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
            families[lineSplit[0]] = Family(lineSplit[0])
        sample = Sample(families[lineSplit[0]], lineSplit[1], lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
    return families
