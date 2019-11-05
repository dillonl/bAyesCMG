import sys
import gzip
sys.path.append('externals/CharGer')
from charger import charger
# import imp
reload(charger)

class VCF:

    def __init__(self, families, vcfFilePath):
        self.vcfFilePath = vcfFilePath
        self.families = families
        self.sampleIdxs = {}
        print('1')
        self.charger = charger.charger()
        # [ vepDone , preVEP , exacDone , clinvarDone ] = self.charger.getInputData( vcf=vcfFilePath )
        # CharGer.getExternalData( clinvar=True , exac=True , vep=True )
        print('hello world')
        '''
        vcf=vcfFile
        self.charger.PVS1( )
        self.charger.PS1( )
        self.charger.PS2( )
        self.charger.PS3( )
        self.charger.PS4( )
        self.charger.PM1( recurrenceThreshold , hotspot3d=clustersFile )
        self.charger.PM2( rareAF )
        self.charger.PM3( )
        self.charger.PM4( )
        self.charger.PM5( )
        self.charger.PM6( )
        self.charger.PP1( )
        self.charger.PP2( )
        self.charger.PP3( minimumEvidence )
        self.charger.PP4( )
        self.charger.PP5( )

        self.charger.BA1( commonAF )
        self.charger.BS1( )
        self.charger.BS2( )
        self.charger.BS3( )
        self.charger.BS4( )
        self.charger.BP1( )
        self.charger.BP2( )
        self.charger.BP3( )
        self.charger.BP4( minimumEvidence )
        self.charger.BP5( )
        self.charger.BP6( )
        self.charger.BP7( )

        self.charger.PSC1( )
        self.charger.PMC1( )
        self.charger.PPC1( )
        self.charger.PPC2( )

        self.charger.BSC1( )
        self.charger.BMC1( )

        print( str( rareAF ) + " < " + str( commonAF ) )
        t4 = time.time()

        self.charger.classify( system="ACMG", scoresMap = scoresMap )
        '''

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
