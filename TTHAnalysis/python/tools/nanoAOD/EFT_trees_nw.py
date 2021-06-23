from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.physicsobjects import _btagWPs as HiggsRecoTTHbtagwps

import ROOT, itertools
from ROOT import *
import numpy as np
import math
import os
import tensorflow as tf

import h5py

from keras import optimizers
from keras.models import load_model

class EFTtrees_nw(Module):
    def __init__(self,label="_Recl", variations=[], cut_BDT_rTT_score = 0.0, btagDeepCSVveto = 'M', doSystJEC=True):

        self.label = label
        self.cut_BDT_rTT_score = cut_BDT_rTT_score
        self.btagDeepCSVveto = btagDeepCSVveto
        self.branches = []
        self.systsJEC = {0:"", 1:"_jesTotalCorrUp", -1:"_jesTotalCorrDown", 2:"_jesTotalUnCorrUp", -2:"_jesTotalUnCorrDown"} if doSystJEC else {0:""}
        if(doSystJEC is True):
            if len(variations):
                self.systsJEC = {0:""}
                for i,var in enumerate(variations):
                    self.systsJEC[i+1]   ="_%sUp"%var
                    self.systsJEC[-(i+1)]="_%sDown"%var
        else: self.systsJEC = {0:""}
        self.nlep = 2
        self.njet = 5
        self.ngenjet = 8


    def setDefault(self, event, jesLabel):
        for jesLabel in self.systsJEC.values():

            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                self.out.fillBranch('%sHadTop%s%s'%(self.label,jesLabel,suffix), -99.)

            # Jets
            for suffix in ["_pt", "_eta", "_phi", "_mass", "_btagdiscr", "_ishadtop"]:
                for iJet in range(self.njet):
                    self.out.fillBranch('%sJet%s%s%s'%(self.label,iJet,jesLabel,suffix)   , -99.)

            self.out.fillBranch('%sTopScore%s'%(self.label,jesLabel)      , -99.)
            self.out.fillBranch('%smet%s'%(self.label,jesLabel)           , -99.)
            self.out.fillBranch('%smet_phi%s'%(self.label,jesLabel)       , -99.)
            
            self.out.fillBranch('%sJet_Higgs_score%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sreco_score_higgs%s'%(self.label,jesLabel), -99.)



            for var in ['DeltaRl0l1',]:
                self.out.fillBranch('%s%s%s'%(self.label,var,jesLabel), -99.)


            for iLep in range(self.nlep):
                self.out.fillBranch('%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , -99.)
                self.out.fillBranch('%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , -99.)
                self.out.fillBranch('%sClosestJetToLep_pt%s%s'%(self.label,iLep,jesLabel) , -99.)
                self.out.fillBranch('%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel) , -99.)
                self.out.fillBranch('%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel) , -99.)
                self.out.fillBranch('%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel) , -99.)


            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                for iLep in range(self.nlep):
                    self.out.fillBranch('%sLep%s%s%s'%(self.label,iLep,jesLabel,suffix)   , -99.)

            self.out.fillBranch('%sMore5_Jets_pt%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sMore5_Jets_eta%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sMore5_Jets_phi%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sMore5_Jets_mass%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%sAll5_Jets_pt%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sAll5_Jets_eta%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sAll5_Jets_phi%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sAll5_Jets_mass%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%sJets_plus_Lep_pt%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sJets_plus_Lep_eta%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sJets_plus_Lep_phi%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%sJets_plus_Lep_mass%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%sdnn_prediction%s'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%s_weight_ctp1%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctp2%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctp3%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctp5%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctp10%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%s_weight_cpt1%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cpt2%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cpt3%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cpt5%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cpt10%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%s_weight_cptb1%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cptb2%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cptb3%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cptb5%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_cptb10%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%s_weight_ctG1%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctG2%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctG3%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctG5%s'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%s_weight_ctG10%s'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG1'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_ctp1'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cpt1'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb1'%(self.label, jesLabel), -99.)


            self.out.fillBranch('%sdnn_prediction%s_weight_ctG2'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_ctp2'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cpt2'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb2'%(self.label, jesLabel), -99.)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG3'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_ctp3'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cpt3'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb3'%(self.label, jesLabel), -99.)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG5'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_ctp5'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cpt5'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb5'%(self.label, jesLabel), -99.)


            self.out.fillBranch('%sdnn_prediction%s_weight_ctG10'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_ctp10'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cpt10'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb10'%(self.label, jesLabel), -99.)

            self.out.fillBranch('%sdnn_prediction%s_weight_cptSM'%(self.label, jesLabel), -99.)



            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG1'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp1'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt1'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb1'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG2'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp2'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt2'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb2'%(self.label,jesLabel), -99.)


            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG3'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp3'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt3'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb3'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG5'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp5'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt5'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb5'%(self.label,jesLabel), -99.)



            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG10'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp10'%(self.label, jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt10'%(self.label,jesLabel), -99.)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb10'%(self.label,jesLabel), -99.)

            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptSM'%(self.label, jesLabel), -99.)

#            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt2'%(self.label,jesLabel), -99.)
#            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt3'%(self.label,jesLabel), -99.)
#            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptSM'%(self.label,jesLabel), -99.)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #model_dnn = load_model('dnn_tagger_new_dr.h5') 
        self.out = wrappedOutputTree

        # Somehow dependent on JES

        for jesLabel in self.systsJEC.values():

           # Leptons and the precomputed hadronic top
            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                for iLep in range(self.nlep):
                    self.out.branch('%sLep%s%s%s'%(self.label,iLep,jesLabel,suffix)   , 'F')
                self.out.branch('%sHadTop%s%s'%(self.label,jesLabel,suffix), 'F')

            # Jets
            for suffix in ["_pt", "_eta", "_phi", "_mass", "_btagdiscr", "_ishadtop"]:
                for iJet in range(self.njet):
                    self.out.branch('%sJet%s%s%s'%(self.label,iJet,jesLabel,suffix)   , 'F')


            for iLep in range(self.nlep):
                self.out.branch('%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_pt%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel) , 'F')

            self.out.branch('%sTopScore%s'%(self.label,jesLabel)      , 'F')
            self.out.branch('%smet%s'%(self.label,jesLabel)           , 'F')
            self.out.branch('%smet_phi%s'%(self.label,jesLabel)       , 'F')

            self.out.branch('%sJet_Higgs_score%s'%(self.label,jesLabel), 'F')
            self.out.branch('%sreco_score_higgs%s'%(self.label,jesLabel), 'F')


            for var in ['DeltaRl0l1',]:
                self.out.branch('%s%s%s'%(self.label,var,jesLabel), 'F')

            self.out.branch('%sMore5_Jets_pt%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sMore5_Jets_eta%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sMore5_Jets_phi%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sMore5_Jets_mass%s'%(self.label,jesLabel), 'F'   )

            self.out.branch('%sAll5_Jets_pt%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sAll5_Jets_eta%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sAll5_Jets_phi%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sAll5_Jets_mass%s'%(self.label,jesLabel), 'F'   )

            self.out.branch('%sJets_plus_Lep_pt%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sJets_plus_Lep_eta%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sJets_plus_Lep_phi%s'%(self.label,jesLabel), 'F'   )
            self.out.branch('%sJets_plus_Lep_mass%s'%(self.label,jesLabel), 'F'   )

            self.out.branch('%sdnn_prediction%s'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s'%(self.label,jesLabel), 'F')

            self.out.branch('%s_weight_ctp1%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctp2%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctp3%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctp5%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctp10%s'%(self.label,jesLabel), 'F')

            self.out.branch('%s_weight_cpt1%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cpt2%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cpt3%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cpt5%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cpt10%s'%(self.label,jesLabel), 'F')

            self.out.branch('%s_weight_cptb1%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cptb2%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cptb3%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cptb5%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_cptb10%s'%(self.label,jesLabel), 'F')

            self.out.branch('%s_weight_ctG1%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctG2%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctG3%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctG5%s'%(self.label,jesLabel), 'F')
            self.out.branch('%s_weight_ctG10%s'%(self.label,jesLabel), 'F')


            self.out.branch('%sdnn_prediction%s_weight_ctG1'%(self.label, jesLabel), 'F')        
            self.out.branch('%sdnn_prediction%s_weight_ctp1'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cptb1'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cpt1'%(self.label, jesLabel), 'F')

            self.out.branch('%sdnn_prediction%s_weight_ctG2'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_ctp2'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cptb2'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cpt2'%(self.label, jesLabel), 'F')


            self.out.branch('%sdnn_prediction%s_weight_ctG3'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_ctp3'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cptb3'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cpt3'%(self.label, jesLabel), 'F')

            self.out.branch('%sdnn_prediction%s_weight_ctG5'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_ctp5'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cptb5'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cpt5'%(self.label, jesLabel), 'F')

            self.out.branch('%sdnn_prediction%s_weight_ctG10'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_ctp10'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cptb10'%(self.label, jesLabel), 'F')
            self.out.branch('%sdnn_prediction%s_weight_cpt10'%(self.label, jesLabel), 'F')

            self.out.branch('%sdnn_prediction%s_weight_cptSM'%(self.label, jesLabel), 'F')


            self.out.branch('%shiggs_reco_mass%s_weight_ctG1'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_ctp1'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cpt1'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cptb1'%(self.label, jesLabel), 'F')

            self.out.branch('%shiggs_reco_mass%s_weight_ctG2'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_ctp2'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cpt2'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cptb2'%(self.label, jesLabel), 'F')

            self.out.branch('%shiggs_reco_mass%s_weight_ctG3'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_ctp3'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cpt3'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cptb3'%(self.label, jesLabel), 'F')

            self.out.branch('%shiggs_reco_mass%s_weight_ctG5'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_ctp5'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cpt5'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cptb5'%(self.label, jesLabel), 'F')

            self.out.branch('%shiggs_reco_mass%s_weight_ctG10'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_ctp10'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cpt10'%(self.label, jesLabel), 'F')
            self.out.branch('%shiggs_reco_mass%s_weight_cptb10'%(self.label, jesLabel), 'F')

            self.out.branch('%shiggs_reco_mass%s_weight_cptSM'%(self.label,jesLabel), 'F')

            
    def analyze(self, event):

        # Some useful input parameters
        year=getattr(event,'year')
        #btagvetoval=HiggsRecoTTHbtagwps['DeepFlav_%d_%s'%(year,self.btagDeepCSVveto)][1]

        nAllLeps = event.nLepGood
        nRecleanedLeps = event.nLepFO_Recl
        recleanedLepsIdxs = event.iLepFO_Recl
        allLeps = Collection(event,"LepGood","nLepGood")
        leps = [allLeps[recleanedLepsIdxs[i]] for i in xrange(nRecleanedLeps)]
        alljets = [x for x in Collection(event,"JetSel_Recl","nJetSel_Recl")]
        
        #allLeps = Collection(event,"LepGood","nLepGood")
        #leps = [allLeps[recleanedLepsIdxs[i]] for i in xrange(nRecleanedLeps)]
        #alljets = [x for x in Collection(event,"JetSel_Recl","nJetSel_Recl")]

        for jesLabel in self.systsJEC.values():

            if len(leps) < 2                      :
                self.setDefault(event, jesLabel)
                continue

            #if leps[0].pt < 20 or leps[1].pt < 10 :
            #    self.setDefault(event, jesLabel)
            #    continue
                
            #if event.nLepTight_Recl > 2           :
            #    self.setDefault(event, jesLabel)
            #    continue
                
            #if leps[0].pdgId*leps[1].pdgId < 0    :
            #    self.setDefault(event, jesLabel)
            #    continue

            #if abs(event.mZ1_Recl-91.2)<10        :
            #    self.setDefault(event, jesLabel)
            #    continue

            #if event.nTauSel_Recl_Tight > 0                           :
            #    self.setDefault(event, jesLabel)
            #    continue
            
            #if not( getattr(event,"nJet25%s_Recl"%jesLabel)>=3 and getattr(event,"nBJetLoose25%s_Recl"%jesLabel)>= 1 and getattr(event,"nBJetMedium25%s_Recl"%jesLabel)>= 0 ):
            #if not(getattr(event,"nJet25%s_Recl"%jesLabel)>=4 and getattr(event,"nBJetLoose25%s_Recl"%jesLabel)>= 1 and getattr(event,"nBJetMedium25%s_Recl"%jesLabel)>= 0):
            #    self.setDefault(event, jesLabel)
            #    continue

            for iLep in range(self.nlep):
                part = leps[iLep].p4()
                #                self.out.fillBranch('%sLep%s%s_pt'%(self.label,iLep,jesLabel), getattr(event, '%sLep%s%s_pt'%(self.label,iLep,jesLabel)) )#part.Pt()  )
                #                self.out.fillBranch('%sLep%s%s_eta'%(self.label,iLep,jesLabel), getattr(event, '%sLep%s%s_eta'%(self.label,iLep,jesLabel)))#part.Eta() )
                #                self.out.fillBranch('%sLep%s%s_phi'%(self.label,iLep,jesLabel), getattr(event, '%sLep%s%s_phi'%(self.label,iLep,jesLabel)))#part.Phi() )
                #                self.out.fillBranch('%sLep%s%s_mass'%(self.label,iLep,jesLabel), getattr(event, '%sLep%s%s_mass'%(self.label,iLep,jesLabel)))#part.M()   )
                self.out.fillBranch('%sLep%s%s_pt'%(self.label,iLep,jesLabel), part.Pt() )
                self.out.fillBranch('%sLep%s%s_eta'%(self.label,iLep,jesLabel), part.Eta() )
                self.out.fillBranch('%sLep%s%s_phi'%(self.label,iLep,jesLabel), part.Phi() )
                self.out.fillBranch('%sLep%s%s_mass'%(self.label,iLep,jesLabel), part.M() )

                self.out.fillBranch('%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel)))
                self.out.fillBranch('%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel)))
                self.out.fillBranch('%sClosestJetToLep_pt%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel)))
                self.out.fillBranch('%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel)))
                self.out.fillBranch('%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel)))
                self.out.fillBranch('%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel) , getattr(event, '%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel)))


            dnn_score = getattr(event,'Hreco_dnn_prediction%s'%(jesLabel))
            mass = getattr(event,'Hreco_higgs_reco_mass%s'%(jesLabel))
            
            self.out.fillBranch('%sDeltaRl0l1%s' %(self.label,jesLabel), getattr(event, '%sDeltaRl0l1%s' %(self.label,jesLabel)))



            self.out.fillBranch('%sMore5_Jets_pt%s'%(self.label,jesLabel), getattr(event, '%sMore5_Jets_pt%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sMore5_Jets_eta%s'%(self.label,jesLabel), getattr(event, '%sMore5_Jets_eta%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sMore5_Jets_phi%s'%(self.label,jesLabel), getattr(event, '%sMore5_Jets_phi%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sMore5_Jets_mass%s'%(self.label,jesLabel), getattr(event, '%sMore5_Jets_mass%s'%(self.label,jesLabel)))

            self.out.fillBranch('%sAll5_Jets_pt%s'%(self.label,jesLabel), getattr(event, '%sAll5_Jets_pt%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sAll5_Jets_eta%s'%(self.label,jesLabel), getattr(event, '%sAll5_Jets_eta%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sAll5_Jets_phi%s'%(self.label,jesLabel), getattr(event, '%sAll5_Jets_phi%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sAll5_Jets_mass%s'%(self.label,jesLabel), getattr(event, '%sAll5_Jets_mass%s'%(self.label,jesLabel)))


            self.out.fillBranch('%sJets_plus_Lep_pt%s'%(self.label,jesLabel), getattr(event, '%sJets_plus_Lep_pt%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sJets_plus_Lep_eta%s'%(self.label,jesLabel), getattr(event, '%sJets_plus_Lep_eta%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sJets_plus_Lep_phi%s'%(self.label,jesLabel), getattr(event, '%sJets_plus_Lep_phi%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sJets_plus_Lep_mass%s'%(self.label,jesLabel), getattr(event, '%sJets_plus_Lep_mass%s'%(self.label,jesLabel)))

            for suffix in ["_pt", "_eta", "_phi", "_mass", "_btagdiscr", "_ishadtop"]:
                for iJet in range(self.njet):
                    self.out.fillBranch('%sJet%s%s%s'%(self.label,iJet,jesLabel,suffix)   , getattr(event, '%sJet%s%s%s'%(self.label,iJet,jesLabel,suffix)))

            self.out.fillBranch('%sTopScore%s'%(self.label,jesLabel)      , getattr(event, '%sTopScore%s'%(self.label,jesLabel)))
            self.out.fillBranch('%smet%s'%(self.label,jesLabel)           , getattr(event, '%smet%s'%(self.label,jesLabel)))
            self.out.fillBranch('%smet_phi%s'%(self.label,jesLabel)       , getattr(event, '%smet_phi%s'%(self.label,jesLabel)))

            self.out.fillBranch('%sJet_Higgs_score%s'%(self.label,jesLabel), getattr(event, '%sJet_Higgs_score%s'%(self.label,jesLabel)))
            self.out.fillBranch('%sreco_score_higgs%s'%(self.label,jesLabel), getattr(event, '%sreco_score_higgs%s'%(self.label,jesLabel)))


            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                self.out.fillBranch('%sHadTop%s%s'%(self.label,jesLabel,suffix), getattr(event, '%sHadTop%s%s'%(self.label,jesLabel,suffix)))




            if(dnn_score == -99 or mass == -99.):
                self.setDefault(event, jesLabel)
                continue

            #cpt reweight
            ctp_1 = getattr(event, 'weight_ttH_EFT1_ctp.dat')            
            cpt_1 = getattr(event, 'weight_ttH_EFT1_cpt.dat')
            ctG_1 = getattr(event, 'weight_ttH_EFT1_ctG.dat')
            cptb_1 = getattr(event, 'weight_ttH_EFT1_cptb.dat')

            ctp_2 = getattr(event, 'weight_ttH_EFT2a_ctp.dat')
            cpt_2 = getattr(event, 'weight_ttH_EFT2a_cpt.dat')
            ctG_2 = getattr(event, 'weight_ttH_EFT2a_ctG.dat')
            cptb_2 = getattr(event, 'weight_ttH_EFT2a_cptb.dat')

            ctp_3 = getattr(event, 'weight_ttH_EFT3a_ctp.dat')
            cpt_3 = getattr(event, 'weight_ttH_EFT3a_cpt.dat')
            ctG_3 = getattr(event, 'weight_ttH_EFT3a_ctG.dat')
            cptb_3 = getattr(event, 'weight_ttH_EFT3a_cptb.dat')

            ctp_5 = getattr(event, 'weight_ttH_EFT5a_ctp.dat')
            cpt_5 = getattr(event, 'weight_ttH_EFT5a_cpt.dat')
            ctG_5 = getattr(event, 'weight_ttH_EFT5a_ctG.dat')
            cptb_5 = getattr(event, 'weight_ttH_EFT5a_cptb.dat')

            ctp_10 = getattr(event, 'weight_ttH_EFT10a_ctp.dat')
            cpt_10 = getattr(event, 'weight_ttH_EFT10a_cpt.dat')
            ctG_10 = getattr(event, 'weight_ttH_EFT10a_ctG.dat')
            cptb_10 = getattr(event, 'weight_ttH_EFT10a_cptb.dat')

#            cpt_1 = getattr(event, 'weight_ttH_EFT1_cpt.dat')#*0.5071
#            cpt_2 = getattr(event, 'weight_ttH_EFT2a_cpt.dat')#*0.5071
#            cpt_3 = getattr(event, 'weight_ttH_EFT3a_cpt.dat')#*0.5071
            cpt_SM = getattr(event, 'weight_sm.dat')#0.5071
       
            self.out.fillBranch('%sdnn_prediction%s'%(self.label, jesLabel), dnn_score)
            self.out.fillBranch('%shiggs_reco_mass%s'%(self.label,jesLabel), mass)

            self.out.fillBranch('%s_weight_ctp1%s'%(self.label,jesLabel), ctp_1)
            self.out.fillBranch('%s_weight_ctp2%s'%(self.label,jesLabel), ctp_2)
            self.out.fillBranch('%s_weight_ctp3%s'%(self.label,jesLabel), ctp_3)
            self.out.fillBranch('%s_weight_ctp5%s'%(self.label,jesLabel), ctp_5)
            self.out.fillBranch('%s_weight_ctp10%s'%(self.label,jesLabel), ctp_10)

            self.out.fillBranch('%s_weight_cpt1%s'%(self.label,jesLabel), cpt_1)
            self.out.fillBranch('%s_weight_cpt2%s'%(self.label,jesLabel), cpt_2)
            self.out.fillBranch('%s_weight_cpt3%s'%(self.label,jesLabel), cpt_3)
            self.out.fillBranch('%s_weight_cpt5%s'%(self.label,jesLabel), cpt_5)
            self.out.fillBranch('%s_weight_cpt10%s'%(self.label,jesLabel), cpt_10)

            self.out.fillBranch('%s_weight_cptb1%s'%(self.label,jesLabel), cptb_1)
            self.out.fillBranch('%s_weight_cptb2%s'%(self.label,jesLabel), cptb_2)
            self.out.fillBranch('%s_weight_cptb3%s'%(self.label,jesLabel), cptb_3)
            self.out.fillBranch('%s_weight_cptb5%s'%(self.label,jesLabel), cptb_5)
            self.out.fillBranch('%s_weight_cptb10%s'%(self.label,jesLabel), cptb_10)

            self.out.fillBranch('%s_weight_ctG1%s'%(self.label,jesLabel), ctG_1)
            self.out.fillBranch('%s_weight_ctG2%s'%(self.label,jesLabel), ctG_2)
            self.out.fillBranch('%s_weight_ctG3%s'%(self.label,jesLabel), ctG_3)
            self.out.fillBranch('%s_weight_ctG5%s'%(self.label,jesLabel), ctG_5)
            self.out.fillBranch('%s_weight_ctG10%s'%(self.label,jesLabel), ctG_10)


            self.out.fillBranch('%sdnn_prediction%s_weight_ctp1'%(self.label, jesLabel), dnn_score*ctp_1)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp1'%(self.label,jesLabel), mass*ctp_1)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG1'%(self.label, jesLabel), dnn_score*ctG_1)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG1'%(self.label,jesLabel), mass*ctG_1)

            self.out.fillBranch('%sdnn_prediction%s_weight_cpt1'%(self.label, jesLabel), dnn_score*cpt_1)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt1'%(self.label,jesLabel), mass*cpt_1)

 
            self.out.fillBranch('%sdnn_prediction%s_weight_cptb1'%(self.label, jesLabel), dnn_score*cptb_1)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb1'%(self.label,jesLabel), mass*cptb_1)
           

#            self.out.fillBranch('%sdnn_prediction%s_weight_cpt2'%(self.label, jesLabel), dnn_score*cpt_2)
#            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt2'%(self.label,jesLabel), mass*cpt_2)


            self.out.fillBranch('%sdnn_prediction%s_weight_cpt2'%(self.label, jesLabel), dnn_score*cpt_2)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt2'%(self.label,jesLabel), mass*cpt_2)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG2'%(self.label, jesLabel), dnn_score*ctG_2)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG2'%(self.label,jesLabel), mass*ctG_2)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctp2'%(self.label, jesLabel), dnn_score*ctp_2)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp2'%(self.label,jesLabel), mass*ctp_2)

            self.out.fillBranch('%sdnn_prediction%s_weight_cptb2'%(self.label, jesLabel), dnn_score*cptb_2)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb2'%(self.label,jesLabel), mass*cptb_2)




            self.out.fillBranch('%sdnn_prediction%s_weight_cpt3'%(self.label, jesLabel), dnn_score*cpt_3)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt3'%(self.label,jesLabel), mass*cpt_3)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG3'%(self.label, jesLabel), dnn_score*ctG_3)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG3'%(self.label,jesLabel), mass*ctG_3)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctp3'%(self.label, jesLabel), dnn_score*ctp_3)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp3'%(self.label,jesLabel), mass*ctp_3)

            self.out.fillBranch('%sdnn_prediction%s_weight_cptb3'%(self.label, jesLabel), dnn_score*cptb_3)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb3'%(self.label,jesLabel), mass*cptb_3)


            self.out.fillBranch('%sdnn_prediction%s_weight_cpt5'%(self.label, jesLabel), dnn_score*cpt_5)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt5'%(self.label,jesLabel), mass*cpt_5)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG5'%(self.label, jesLabel), dnn_score*ctG_5)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG5'%(self.label,jesLabel), mass*ctG_5)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctp5'%(self.label, jesLabel), dnn_score*ctp_5)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp5'%(self.label,jesLabel), mass*ctp_5)

            self.out.fillBranch('%sdnn_prediction%s_weight_cptb5'%(self.label, jesLabel), dnn_score*cptb_5)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb5'%(self.label,jesLabel), mass*cptb_5)


            self.out.fillBranch('%sdnn_prediction%s_weight_cpt10'%(self.label, jesLabel), dnn_score*cpt_10)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cpt10'%(self.label,jesLabel), mass*cpt_10)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctG10'%(self.label, jesLabel), dnn_score*ctG_10)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctG10'%(self.label,jesLabel), mass*ctG_10)

            self.out.fillBranch('%sdnn_prediction%s_weight_ctp10'%(self.label, jesLabel), dnn_score*ctp_10)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_ctp10'%(self.label,jesLabel), mass*ctp_10)


            self.out.fillBranch('%sdnn_prediction%s_weight_cptb10'%(self.label, jesLabel), dnn_score*cptb_10)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptb10'%(self.label,jesLabel), mass*cptb_10)

            self.out.fillBranch('%sdnn_prediction%s_weight_cptSM'%(self.label, jesLabel), dnn_score*cpt_SM)
            self.out.fillBranch('%shiggs_reco_mass%s_weight_cptSM'%(self.label,jesLabel), mass*cpt_SM)



        return True

eft_trees_nw = lambda : EFTtrees_nw(label='Hreco_')
#ttH_2lss_dnn_pt_regression = lambda : Class_ttH_2lss_dnn_pt_regression(label='Hreco_',
#                                                         btagDeepCSVveto = 'M')
