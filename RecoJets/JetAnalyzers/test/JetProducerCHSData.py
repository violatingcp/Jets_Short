# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register ('useData',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.int,
		  "Run this on real data")

options.register ('forceCheckClosestZVertex',
		  False,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.int,
		  "Force the check of the closest z vertex")

if not options.useData :
	inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
else :
	inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
	
# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff.
#from PhysicsTools.PatAlgos.tools.pfTools import *
#postfix = "PFlow"
#usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfix,
#	  jetCorrections=inputJetCorrLabel, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
#if not options.forceCheckClosestZVertex :
#        process.pfPileUpPFlow.checkClosestZVertex = False

process.load('CommonTools.ParticleFlow.pfNoPileUp_cff')
process.load('CommonTools.ParticleFlow.pfParticleSelection_cff')
process.load('CommonTools.ParticleFlow.pfPhotons_cff')
process.load('CommonTools.ParticleFlow.pfElectrons_cff')
process.load('CommonTools.ParticleFlow.pfMuons_cff')
process.load('CommonTools.ParticleFlow.pfJets_cff')
process.load('CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi')
process.load('CommonTools.ParticleFlow.TopProjectors.pfNoElectron_cfi')
process.load('CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi')
process.pfPileUp.PFCandidates          = 'particleFlow'
process.pfNoPileUp.bottomCollection    = 'particleFlow'
process.pfPileUpIso.PFCandidates       = 'particleFlow'
process.pfNoPileUpIso.bottomCollection = 'particleFlow'
process.pfPileUp.Enable                = True
process.pfPileUp.Vertices              = 'goodOfflinePrimaryVertices'
process.pfPileUp.checkClosestZVertex   = cms.bool(False)
process.pfJets.doAreaFastjet           = True
process.pfJets.doRhoFastjet            = False

process.ak5PFJetsCHS = cms.EDProducer("FastjetJetProducer",
				         Active_Area_Repeats = cms.int32(1),
				         doAreaFastjet = cms.bool(True),
				         voronoiRfact = cms.double(-0.9),
				         maxBadHcalCells = cms.uint32(9999999),
				         doAreaDiskApprox = cms.bool(False),
				         maxRecoveredEcalCells = cms.uint32(9999999),
				         jetType = cms.string('PFJet'),
				         minSeed = cms.uint32(14327),
				         Ghost_EtaMax = cms.double(5.0),
				         doRhoFastjet = cms.bool(False),
				         jetAlgorithm = cms.string('AntiKt'),
				         nSigmaPU = cms.double(1.0),
				         GhostArea = cms.double(0.01),
				         Rho_EtaMax = cms.double(4.4),
				         maxBadEcalCells = cms.uint32(9999999),
				         useDeterministicSeed = cms.bool(True),
				         doPVCorrection = cms.bool(False),
				         maxRecoveredHcalCells = cms.uint32(9999999),
				         rParam = cms.double(0.5),
				         maxProblematicHcalCells = cms.uint32(9999999),
				         doOutputJets = cms.bool(True),
				         src = cms.InputTag("pfNoElectron"),
				         inputEtMin = cms.double(0.0),
				         srcPVs = cms.InputTag("goodOfflinePrimaryVertices"),
				         jetPtMin = cms.double(3.0),
				         radiusPU = cms.double(0.5),
				         maxProblematicEcalCells = cms.uint32(9999999),
				         doPUOffsetCorr = cms.bool(False),
				         inputEMin = cms.double(0.0)
				     )


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
         	"PrimaryVertexObjectFilter",
	        filterParams = pvSelector.clone( maxZ = cms.double(24.0),
						 minNdof = cms.double(4.0) # this is >= 4
						 ),
	        src=cms.InputTag("offlinePrimaryVertices")
	        )


process.load("RecoJets.Configuration.GenJetParticles_cff")

process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
					         src = cms.InputTag("genParticles"),
					         ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012,
										         2000014, 2000016, 1000039, 5100039, 4000012,
										         4000014, 4000016, 9900012, 9900014, 9900016,
										         39),
					         partonicFinalState = cms.bool(False),
					         excludeResonances = cms.bool(True),
					         excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
					         tausAsJets = cms.bool(False)
					     )


from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak5GenJets     = ak5GenJets.clone()
process.ak5GenJetsNoNu = ak5GenJets.clone(src = cms.InputTag("genParticlesForJetsNoNu"))
process.ca8GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))
process.ak8GenJetsNoNu = ak5GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                            src = cms.InputTag("genParticles"),
                                            select = cms.vstring(
    "drop  *"
    ,"keep status = 3" #keeps  particles from the hard matrix element
    ,"keep (abs(pdgId) >= 11 & abs(pdgId) <= 16) & status = 1" #keeps e/mu and nus with status 1
    ,"keep (abs(pdgId)  = 15) & status = 3" #keeps taus
    )
                                            )

from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *

from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJetsPFlow = ca4PFJets.clone(
        rParam = cms.double(0.8),
            src = cms.InputTag('particleFlow'),
            doAreaFastjet = cms.bool(True),
            doRhoFastjet = cms.bool(True),
            Rho_EtaMax = cms.double(6.0),
            Ghost_EtaMax = cms.double(7.0)
            )

process.ca8PFJetsPFlowCHS = ca4PFJets.clone(
        rParam = cms.double(0.8),
            src = cms.InputTag('pfNoElectron'),
            doAreaFastjet = cms.bool(True),
            doRhoFastjet = cms.bool(True),
            Rho_EtaMax = cms.double(6.0),
            Ghost_EtaMax = cms.double(7.0)
            )

from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.ak5TrimmedPFlow = ak5PFJetsTrimmed.clone(
        src = process.ak5PFJetsCHS.src,
            doAreaFastjet = cms.bool(True)
            )

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.ak5FilteredPFlow = ak5PFJetsFiltered.clone(
        src = process.ak5PFJetsCHS.src,
            doAreaFastjet = cms.bool(True)
            )

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ak5PrunedPFlow = ak5PFJetsPruned.clone(
        src = process.ak5PFJetsCHS.src,
            doAreaFastjet = cms.bool(True)
            )
###############################
###### AK 0.8 jets groomed ####
###############################

process.ak8TrimmedPFlow = process.ak5TrimmedPFlow.clone(
            src = process.ak5PFJetsCHS.src,
                    rParam = cms.double(0.8)
                )

process.ak8FilteredPFlow = process.ak5FilteredPFlow.clone(
            src = process.ak5PFJetsCHS.src,
                    rParam = cms.double(0.8)
                    )

process.ak8PrunedPFlow = process.ak5PrunedPFlow.clone(
            src = process.ak5PFJetsCHS.src,
                    rParam = cms.double(0.8)
                )
###############################
###### CA8 Pruning Setup ######
###############################


# Pruned PF Jets
process.caPrunedPFlow = process.ak5PrunedPFlow.clone(
            jetAlgorithm = cms.string("CambridgeAachen"),
	    rParam       = cms.double(0.8)
            )


process.caPrunedGen = process.ca8GenJetsNoNu.clone(
            SubJetParameters,
                    usePruning = cms.bool(True),
                    useExplicitGhosts = cms.bool(True),
                    writeCompound = cms.bool(True),
                    jetCollInstanceName=cms.string("SubJets")
            )

###############################
###### CA8 Filtered Setup #####
###############################


# Filtered PF Jets
process.caFilteredPFlow = ak5PFJetsFiltered.clone(
            src = cms.InputTag('particleFlow'),
	    jetAlgorithm = cms.string("CambridgeAachen"),
	    rParam       = cms.double(1.2),
	    writeCompound = cms.bool(True),
	    doAreaFastjet = cms.bool(True),
	    jetPtMin = cms.double(100.0)
            )

process.caFilteredPFlowCHS = ak5PFJetsFiltered.clone(
            src = cms.InputTag('pfNoElectron'),
	    jetAlgorithm = cms.string("CambridgeAachen"),
	    rParam       = cms.double(1.2),
	    writeCompound = cms.bool(True),
	    doAreaFastjet = cms.bool(True),
	    jetPtMin = cms.double(100.0)
            )

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsMassDropFiltered
process.caMassDropFilteredPFlow = ak5PFJetsMassDropFiltered.clone(
	     src = cms.InputTag('particleFlow'),
	     jetAlgorithm = cms.string("CambridgeAachen"),
	     rParam       = cms.double(1.2),
	     writeCompound = cms.bool(True),
	     doAreaFastjet = cms.bool(True),
	     jetPtMin = cms.double(100.0)
            )

process.caMassDropFilteredPFlowCHS = ak5PFJetsMassDropFiltered.clone(
	     src = cms.InputTag('pfNoElectron'),
	     jetAlgorithm = cms.string("CambridgeAachen"),
	     rParam       = cms.double(1.2),
	     writeCompound = cms.bool(True),
	     doAreaFastjet = cms.bool(True),
	     jetPtMin = cms.double(100.0)
            )


process.caFilteredGenJetsNoNu = process.ca8GenJetsNoNu.clone(
            nFilt = cms.int32(2),
                    rFilt = cms.double(0.3),
                    useFiltering = cms.bool(True),
                    useExplicitGhosts = cms.bool(True),
                    writeCompound = cms.bool(True),
                    rParam       = cms.double(1.2),
                    jetCollInstanceName=cms.string("SubJets"),
                    jetPtMin = cms.double(100.0)
            )

process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdProducer.jets     = "ak5PFJets"
process.pileupJetIdProducer.vertexes = "goodOfflinePrimaryVertices"

process.load('RecoTauTag.Configuration.RecoPFTauTag_cff')


process.skimmuons =  cms.EDFilter("MuonSelector",
				 src = cms.InputTag('muons'),
				 cut = cms.string(    "(isTrackerMuon) && abs(eta) < 2.5 && pt > 20"+#17. "+
				                      "&& isPFMuon"+
						      "&& globalTrack.isNonnull"+
						      "&& innerTrack.hitPattern.numberOfValidPixelHits > 0"+
						      "&& innerTrack.normalizedChi2 < 10"+
						      "&& numberOfMatches > 0"+
						      "&& innerTrack.hitPattern.numberOfValidTrackerHits>5"+
						      "&& globalTrack.hitPattern.numberOfValidHits>0"+
						      "&& (pfIsolationR03.sumChargedHadronPt+pfIsolationR03.sumNeutralHadronEt+pfIsolationR03.sumPhotonEt)/pt < 0.3"+
						      "&& abs(innerTrack().dxy)<2.0"
						      ),
				 filter = cms.bool(False)
				 )




process.MuonSkim = cms.EDFilter("CandViewCountFilter",
				src = cms.InputTag('skimmuons'),
				minNumber = cms.uint32(2)
				)


process.produce = cms.Sequence(
	    process.skimmuons*
	    process.MuonSkim*
	    #process.scrapingVeto*
            #process.HBHENoiseFilter*
	    process.goodOfflinePrimaryVertices*
	    process.pfNoPileUpSequence*
	    process.pfParticleSelectionSequence*
	    process.pfPhotonSequence*
	    process.pfMuonSequence*
	    process.pfNoMuon*
	    process.pfElectronSequence*
	    process.pfNoElectron*
	    #process.genParticlesForJets*
	    #process.genParticlesForJetsNoNu*
	    #process.ak5GenJets*
	    #process.ak5GenJetsNoNu*
            #process.ca8GenJetsNoNu*
            #process.ak8GenJetsNoNu*
            #process.caFilteredGenJetsNoNu*
	    #process.flavorHistorySeq*
            #process.prunedGenParticles*
            #process.caPrunedGen*
	    process.ak5PFJetsCHS*
	    process.ca8PFJetsPFlow*
	    process.ca8PFJetsPFlowCHS*
	    process.caPrunedPFlow*
	    process.PFTau
	    )

process.p0 = cms.Path(	process.produce	)
 
process.source.fileNames = ['/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/20000/FEF2A4B5-2082-E211-BC3E-E0CB4EA0A904.root']
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.out.SelectEvents.SelectEvents = cms.vstring('p0')
process.out.fileName = cms.untracked.string('output.root')
process.out.dropMetaData = cms.untracked.string("DROPPED")
#process.out.outputCommands   = cms.untracked.vstring("keep *")
process.out.outputCommands = [
	        'keep *_goodOfflinePrimaryVertices*_*_*',
		'keep recoCaloJets_*_*_*',
		'keep recoJPTJets_*_*_*',
		'keep recoPFJets_*_*_*',
		'keep *_*_*_PAT',
#	        'keep *_TriggerResults_*_*',
#	        'keep *_hltTriggerSummaryAOD_*_*',
		'keep *_prunedGenParticles_*_*',
		'drop recoBasicJets_*_*_*',
	        'drop *_*Lite_*_*',
		'keep *_*_rho_*',
		'keep recoElectrons_*_*_*',
	        'keep recoMuons_*_*_*',
		'keep *_offlineBeamSpot_*_*',
#	        'keep *_allConversions_*_*',
		'keep *_pfNoElectron*_*_*',
		'keep recoPFCandidates_particleFlow*_*_*',
		'keep recoTracks_*_*_*',
		'keep *_offlinePrimaryVertices_*_*',
		'keep recoGsfElectrons_*_*_*',
		'keep *_caPrunedPFlow_*_*',
		'keep reco*MET*_*_*_*',
		'keep recoGsfTracks_*_*_*',
		'keep recoGsfElectron*_*_*_*'
	        ]

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

#open('junk.py','w').write(process.dumpPython())
