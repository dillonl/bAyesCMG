'''
countPVst = variant.getPVstCounts()
        countPSt = variant.getPStCounts()
        countPM = variant.getPMCounts()
        countPSu = variant.getPSuCounts()
        oddsPathEquation = self.oddsPathogenicity ** ((countPSu / (self.exponent ** 3)) + (countPM / (self.exponent ** 2)) + (countPSt / (self.exponent ** 1)) + (countPVst / 1))

        countBSt = variant.getBStCounts()
        countBSu = variant.getBSuCounts()
'''

class Variant:

    def __init__(self, variantLine):
        self.line = variantLine

    def getPStCounts(self):
        return 0

    def getPMCounts(self):
        return 0

    def getPSuCounts(self):
        return 0

    def getPVstCounts(self):
        return 0

    def getBStCounts(self):
        return 0

    def getBSuCounts(self):
        return 0
