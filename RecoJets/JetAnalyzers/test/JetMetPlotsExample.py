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
            
#isMC = False
isMC = True

##   ____             __ _                       _     _           
##  / ___|___  _ __  / _(_) __ _ _   _ _ __ __ _| |__ | | ___  ___ 
## | |   / _ \| '_ \| |_| |/ _` | | | | '__/ _` | '_ \| |/ _ \/ __|
## | |__| (_) | | | |  _| | (_| | |_| | | | (_| | |_) | |  __/\__ \
##  \____\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|_.__/|_|\___||___/
##                         |___/                                   

NJetsToKeep = 2
JetPtMin    = 20
UseCHS=''
#UseCHS='CHS'
CaloJetCollection = 'ak5CaloJets'
PFJetCollection   = 'ak5PFJets'+UseCHS
PFJetCollectionCorr   = 'ak5PFJetsCorr'
JPTJetCollection  = 'JetPlusTrackZSPCorJetAntiKt5'
GenJetCollection  = 'ak5GenJets'
PFJetCollectionTight   = "tightPFJetsPFlow"
#CAJetCollection        = "ca8PFJetsPFlow"
CAJetCollection        = "ca8PFJetsPFlow"
CAPrunedJetCollection  = "caPrunedPFlow"
GenJetCollection       = "ak5GenJetsNoNu"
JetCorrection          = "ak5PFL1FastL2L3"

PlotSuffix = "_Data"+UseCHS
inputFile = 'file:output.root'

if isMC:
  PlotSuffix = "_MC"+UseCHS
  inputFile ='file:output.root'

if not isMC:
  JetCorrection += "Residual"


##
## (_)_ __   ___| |_   _  __| | ___  ___ 
## | | '_ \ / __| | | | |/ _` |/ _ \/ __|
## | | | | | (__| | |_| | (_| |  __/\__ \
## |_|_| |_|\___|_|\__,_|\__,_|\___||___/
    
process = cms.Process("Ana")
#############   Format MessageLogger #################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10



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
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('patTuple.root'),
#                               outputCommands = cms.untracked.vstring('keep *')
#                               )
#process.outpath = cms.EndPath(process.out)
##
##
##
process.ak5PFJetsCorr = cms.EDProducer('PFJetCorrectionProducer',
                                       src = cms.InputTag(PFJetCollection),
                                       correctors = cms.vstring(JetCorrection) # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
                                       )


##   ______     _      _      _
##  |__  __|__ | |_   (_) __| |
##  _  | |/ _ \| __/  | |/ _` |
## | \_| |  __/| |_   | | (_| |
##  \____|\___| \__|  |_|\__,_|
##
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.tightPFJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorBasicFilter",
                                         filterParams = pfJetIDSelector.clone(quality=cms.string("TIGHT")),
                                         src = cms.InputTag(PFJetCollectionCorr)
                                         )


process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdProducer.jets     = PFJetCollectionCorr
process.pileupJetIdProducer.vertexes = "goodOfflinePrimaryVertices"
process.pileupJetIdProducer.residualsTxt  = cms.FileInPath("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml")
##  ____  _       _       
## |  _ \| | ___ | |_ ___ 
## | |_) | |/ _ \| __/ __|
## |  __/| | (_) | |_\__ \
## |_|   |_|\___/ \__|___/

#############   Calo Jets  ###########################
process.calo = cms.EDAnalyzer("CaloJetPlotsExample",
    JetAlgorithm  = cms.string(CaloJetCollection),
    HistoFileName = cms.string('CaloJetPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                          
)
#############   PF Jets    ###########################
process.pf = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(PFJetCollectionCorr),
    HistoFileName = cms.string('PFJetPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag("pileupJetIdProducer","fullDiscriminant"),
    PUJetId           = cms.InputTag("pileupJetIdProducer","fullId"),
    JetPtMin          = cms.double(JetPtMin),                          
    PassPUId          = cms.bool(False)                          
)
#############   PF Jets w/PUJet Id   ###########################
process.pfpuid = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(PFJetCollectionCorr),
    HistoFileName = cms.string('PFJetPlotsPUIdExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag("pileupJetIdProducer","fullDiscriminant"),
    PUJetId           = cms.InputTag("pileupJetIdProducer","fullId"),
    JetPtMin          = cms.double(JetPtMin),                          
    PassPUId          = cms.bool(True)                          
)
#############   PF Jets, Tight Jet ID  ################
process.pfTight = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(PFJetCollectionTight),
    HistoFileName = cms.string('PFJetTightPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                          
)

#############   PF Jets, No Corrections    ###########
process.pfUncorr = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(PFJetCollection),
    HistoFileName = cms.string('PFJetUncorrPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    jecLevels     = cms.string("Uncorrected"),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                          
)



#############   JPT Jets    ###########################
process.jpt = cms.EDAnalyzer("JPTJetPlotsExample",
    JetAlgorithm  = cms.string(JPTJetCollection),
    HistoFileName = cms.string('JPTJetPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                          
)
#############   Gen Jets   ###########################
process.gen = cms.EDAnalyzer("GenJetPlotsExample",
     JetAlgorithm  = cms.string(GenJetCollection),
     HistoFileName = cms.string('GenJetPlotsExample'+PlotSuffix+'.root'),
     NJets         = cms.int32(NJetsToKeep),
     PUJetDiscriminant = cms.InputTag(""),
     PUJetId           = cms.InputTag(""),
     JetPtMin          = cms.double(JetPtMin),
     PassPUId          = cms.bool(False)                          
)

#############   Cambridge-Aachen Jets R=0.8 ###########################
process.ca = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(CAJetCollection),
    HistoFileName = cms.string('CAJetPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),                          
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                          
)


#############   Cambridge-Aachen Jets R=0.8, With Pruning ############
process.caPruned = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string(CAPrunedJetCollection),
    HistoFileName = cms.string('CAPrunedJetPlotsExample'+PlotSuffix+'.root'),
    NJets         = cms.int32(NJetsToKeep),
    PUJetDiscriminant = cms.InputTag(""),
    PUJetId           = cms.InputTag(""),
    JetPtMin          = cms.double(JetPtMin),
    PassPUId          = cms.bool(False)                                                            
)

##   _   _        _     
##  | \_/ |  ___ | |_   
##  |     | / _ \| __/  
##  |     ||  __/| |_   
##  |_/\/\| \___|\___|  
##
#Pileup Reduced MET
#Load the conditions (needed for the JEC applied inside the MVA Met)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.GlobalTag.globaltag = 'START53_V22::All'
#Load teh MVA Met
process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
process.load('RecoMET.METPUSubtraction.noPileUpPFMET_cff')

#Specify the leptons
process.pfMEtMVA.srcLeptons =  cms.VInputTag("isomuons")#,"isoelectrons","isotaus")

#Type1 Corrections
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
process.load("JetMETCorrections.Type1MET.correctedMet_cff")
process.corrPfMetType1.jetCorrLabel = cms.string(JetCorrection)

#MET Filters
process.load("RecoMET.METFilters.metFilters_cff")
#Special Set of MET filters that will work with the dataset (THIS IS NOT ALL OF THE FILTERS!!!!!)
#See Here for all https://cmssdt.cern.ch/SDT/lxr/source/RecoMET/METFilters/python/metFilters_cff.py
process.MetFilter = cms.Sequence(process.hcalLaserEventFilter *
                                 process.goodVertices * process.trackingFailureFilter )


#############   caloMET ############
process.caloMet = cms.EDAnalyzer("CaloMETPlotsExample",
    MetAlgorithm  = cms.string("corMetGlobalMuons"),
    HistoFileName = cms.string('CaloMetPlotsExample'+PlotSuffix+'.root'),
    useRecoil     = cms.bool(True),
    recoilLeptons = cms.string("isomuons")                           
)

#############   PFMET ############
process.PFMet = cms.EDAnalyzer("PFMETPlotsExample",
    MetAlgorithm  = cms.string("pfMet"),
    HistoFileName = cms.string('PFMetPlotsExample'+PlotSuffix+'.root'),
    useRecoil     = cms.bool(True),
    recoilLeptons = cms.string("isomuons")                           
)

#############   NoPUMET ############
process.NoPUMet = cms.EDAnalyzer("PFMETPlotsExample",
    MetAlgorithm  = cms.string("pfMEtNoPU"),
    HistoFileName = cms.string('PFMetNoPUPlotsExample'+PlotSuffix+'.root'),
    useRecoil     = cms.bool(True),
    recoilLeptons = cms.string("isomuons")                           
)

#############   MVAMET ############
process.MVAMet = cms.EDAnalyzer("PFMETPlotsExample",
    MetAlgorithm  = cms.string("pfMEtMVA"),
    HistoFileName = cms.string('PFMVAMetPlotsExample'+PlotSuffix+'.root'),
    useRecoil     = cms.bool(True),
    recoilLeptons = cms.string("isomuons")                           
)

#############   PFMet Type 1 correction ############
process.PFMetType1 = cms.EDAnalyzer("PFMETPlotsExample",
    MetAlgorithm  = cms.string("pfMetT1"),
    HistoFileName = cms.string('PFMetT1PlotsExample'+PlotSuffix+'.root'),
    useRecoil     = cms.bool(True),
    recoilLeptons = cms.string("isomuons")                           
)


#############   Path       ###########################
##  ____       _   _     
## |  _ \ __ _| |_| |__  
## | |_) / _` | __| '_ \ 
## |  __/ (_| | |_| | | |
## |_|   \__,_|\__|_| |_|

##process.myseq = cms.Sequence(process.calo*process.pf*process.jpt*process.gen)
process.myseq = cms.Sequence(process.ak5PFJetsCorr* process.pileupJetIdProducer*  process.pfMEtMVAsequence * #process.noPileUpPFMEtSequence *
#                             process.MetFilter *
                             process.pf *process.pfpuid* process.gen  * process.pfUncorr * 
                             process.tightPFJetsPFlow * process.pfTight*  
                             process.ca * process.caPruned * 
#                             process.correctionTermsPfMetType1Type2*process.pfMetT1*process.PFMetType1*
                             process.caloMet*process.PFMet * process.MVAMet #* process.NoPUMet
                             )
  
if not isMC: 
  process.myseq.remove ( process.gen )

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
