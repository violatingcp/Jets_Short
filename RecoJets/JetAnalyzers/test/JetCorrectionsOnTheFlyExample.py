# PYTHON configuration file for class: JetPlotsExample
# Description:  Example of simple EDAnalyzer for jets.
# Author: K. Kousouris
# Date:  25 - August - 2008
# Modified: Kalanand Mishra
# Date:  11 - January - 2011 (for CMS Data Analysis School jet exercise)


import FWCore.ParameterSet.Config as cms

##  ____        _                       __  __  ____ 
## |  _ \  __ _| |_ __ _    ___  _ __  |  \/  |/ ___|
## | | | |/ _` | __/ _` |  / _ \| '__| | |\/| | |    
## | |_| | (_| | || (_| | | (_) | |    | |  | | |___ 
## |____/ \__,_|\__\__,_|  \___/|_|    |_|  |_|\____|
            
isMC = True
##isMC = True

##   ____             __ _                       _     _           
##  / ___|___  _ __  / _(_) __ _ _   _ _ __ __ _| |__ | | ___  ___ 
## | |   / _ \| '_ \| |_| |/ _` | | | | '__/ _` | '_ \| |/ _ \/ __|
## | |__| (_) | | | |  _| | (_| | |_| | | | (_| | |_) | |  __/\__ \
##  \____\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|_.__/|_|\___||___/
##                         |___/

PFJetCollection   = 'ak5PFJets'

PlotSuffix = "_Data"
inputFile = 'file:ttbsm_53x_data_1_1_8sY.root'

if isMC:
  PlotSuffix = "_MC"
  jecLevels = [
    'GR_R_52_V7_L1FastJet_AK5PFchs.txt',
    'GR_R_52_V7_L2Relative_AK5PFchs.txt',
    'GR_R_52_V7_L3Absolute_AK5PFchs.txt'
  ]
else :
  jecLevels = [
    'GR_R_52_V7_L1FastJet_AK5PFchs.txt',
    'GR_R_52_V7_L2Relative_AK5PFchs.txt',
    'GR_R_52_V7_L3Absolute_AK5PFchs.txt',
    'GR_R_52_V7_L2L3Residual_AK5PFchs.txt'
  ]


if isMC:
  PlotSuffix = "_MC"  
  inputFile ='file:ttbsm_52x_mc_128_1_lmR.root'
  

##  _            _           _           
## (_)_ __   ___| |_   _  __| | ___  ___ 
## | | '_ \ / __| | | | |/ _` |/ _ \/ __|
## | | | | | (__| | |_| | (_| |  __/\__ \
## |_|_| |_|\___|_|\__,_|\__,_|\___||___/
    
process = cms.Process("Ana")
#############   Format MessageLogger #################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.TFileService = cms.Service("TFileService",
  fileName = cms.string('jetCorrectionsOnTheFlyExample.root')
)


##  ____             _ ____                           
## |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
## | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
## |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
## |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
                                                   

#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFile)
)
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")


##  ____  _       _       
## |  _ \| | ___ | |_ ___ 
## | |_) | |/ _ \| __/ __|
## |  __/| | (_) | |_\__ \
## |_|   |_|\___/ \__|___/

#############   PF Jets    ###########################
process.pf = cms.EDAnalyzer("JetCorrectionsOnTheFly",
    jetSrc = cms.InputTag(PFJetCollection),
    rhoSrc = cms.InputTag('kt6PFJets', 'rho'),
    pvSrc  = cms.InputTag('goodOfflinePrimaryVertices'),
    jecPayloadNames = cms.vstring( jecLevels ),
    jecUncName = cms.string('GR_R_52_V7_Uncertainty_AK5PFchs.txt') 
)



#############   Path       ###########################
##  ____       _   _     
## |  _ \ __ _| |_| |__  
## | |_) / _` | __| '_ \ 
## |  __/ (_| | |_| | | |
## |_|   \__,_|\__|_| |_|

process.myseq = cms.Sequence(process.pf )
process.p = cms.Path( process.myseq)


if isMC:
  process.source.fileNames = [
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_5.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_501.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_502.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_504.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_508.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_51.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_510.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_511.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_512.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_513.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_514.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_515.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_516.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_518.root",
    "/store/cmst3/user/pharris/cmsdas/DY/PFAOD_52.root"
    ]

if not isMC:
  process.source.fileNames = [
  "/store/cmst3/user/pharris/cmsdas/Data/5E27D95A-3582-E211-A762-001EC9D8A8D0.root",
  "/store/cmst3/user/pharris/cmsdas/Data/5E61832F-2F84-E211-87CA-E0CB4E19F9B8.root",
  "/store/cmst3/user/pharris/cmsdas/Data/5EA2C718-0E82-E211-BE24-00259073E504.root",
  "/store/cmst3/user/pharris/cmsdas/Data/60010BC1-2A82-E211-B3AB-485B39800C3D.root",
  "/store/cmst3/user/pharris/cmsdas/Data/603D96E7-DD81-E211-8ED1-485B3989725C.root",
  "/store/cmst3/user/pharris/cmsdas/Data/6073DE17-BB82-E211-B458-E0CB4E1A1195.root",
  "/store/cmst3/user/pharris/cmsdas/Data/60FD6F66-D882-E211-BC57-00259073E3A8.root",
  "/store/cmst3/user/pharris/cmsdas/Data/64F018F1-8082-E211-AF77-00259074AEAC.root",
  "/store/cmst3/user/pharris/cmsdas/Data/6A7EA788-A882-E211-956F-485B39800C0F.root",
  "/store/cmst3/user/pharris/cmsdas/Data/74631A69-1E84-E211-AED4-002590747DE2.root",
  "/store/cmst3/user/pharris/cmsdas/Data/74D04754-1D82-E211-920B-00261834B527.root",
  "/store/cmst3/user/pharris/cmsdas/Data/7679F804-7596-E211-BC26-002590775158.root",
  "/store/cmst3/user/pharris/cmsdas/Data/767BAF25-9682-E211-B715-001EC9D8BDE7.root",
  "/store/cmst3/user/pharris/cmsdas/Data/78A52856-3E82-E211-ABAC-00259073E44E.root",
  "/store/cmst3/user/pharris/cmsdas/Data/7AB3F26E-7382-E211-B801-00259073E3AC.root",
  "/store/cmst3/user/pharris/cmsdas/Data/82E2D888-2684-E211-ABE8-485B39800BD5.root",
  ]
