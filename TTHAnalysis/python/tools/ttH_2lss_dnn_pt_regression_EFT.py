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

class Class_ttH_2lss_dnn_pt_regression_EFT(Module):
    def __init__(self,label="_Recl", variations=[], cut_BDT_rTT_score = 0.0, btagDeepCSVveto = 'M', doSystJEC=True):
        def loss_MSEDeltaVar(y_true, y_pred):
            y_true = tf.cast(y_true,tf.float32)
            y_pred = tf.cast(y_pred,tf.float32)
            y_true_mean = tf.reduce_mean(y_true)
            y_pred_mean = tf.reduce_mean(y_pred)
            base = tf.reduce_mean((y_true-y_pred)**2)
            var_true = tf.reduce_mean((y_true-y_true_mean)**2)
            var_pred = tf.reduce_mean((y_pred-y_pred_mean)**2)
            var_diff = abs(var_true - var_pred)
            val = base*var_diff
            return val
        #self.genpar = Collection(event,"GenPart","nGenPart")
             
        #print self.genpar
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
        self.model_dnn = load_model(os.path.join(os.environ["CMSSW_BASE"], "src/CMGTools/TTHAnalysis/data/regressionMVA/dnn_tagger_new_dr_real.h5"))
        self.model_regression = load_model(os.path.join(os.environ["CMSSW_BASE"], "src/CMGTools/TTHAnalysis/data/regressionMVA/dnn_trained_2lss.h5"), custom_objects={'loss_MSEDeltaVar': loss_MSEDeltaVar})#os.p\


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #model_dnn = load_model('dnn_tagger_new_dr.h5') 
        self.out = wrappedOutputTree
        self.model_dnn.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

        # Somehow dependent on JES

        for jesLabel in self.systsJEC.values():

            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                for iLep in range(self.nlep):
                    self.out.branch('%sLep%s%s%s'%(self.label,iLep,jesLabel,suffix)   , 'F')          
            # Counters
            self.out.branch('%sdnn_prediction%s'%(self.label, jesLabel)   , 'F')    
            # Gen level, the labels
            self.out.branch('%sHTXS_Higgs_pt'%(self.label) , 'F')
            self.out.branch('%shiggs_reco_mass%s'%(self.label,jesLabel), 'F')

    def setDefault(self, event, jesLabel):

        for suffix in ["_pt", "_eta", "_phi", "_mass"]:
            for iLep in range(self.nlep):
                self.out.fillBranch('%sLep%s%s%s'%(self.label,iLep,jesLabel,suffix)   , -99.)

        self.out.fillBranch('%sdnn_prediction%s'%(self.label, jesLabel)   , -99.)
        # Gen level, the labels
        self.out.fillBranch('%sHTXS_Higgs_pt'%(self.label) , -99.)
            
    def buildHadronicTop(self, event, score, alljets, jesLabel):
        HadTop=None
        if score>self.cut_BDT_rTT_score:
            j1top = int(getattr(event,"BDThttTT_eventReco_iJetSel1%s"%jesLabel))
            j2top = int(getattr(event,"BDThttTT_eventReco_iJetSel2%s"%jesLabel))
            j3top = int(getattr(event,"BDThttTT_eventReco_iJetSel3%s"%jesLabel))
            # Build hadronic top
            top1 = ROOT.TLorentzVector(); top1.SetPtEtaPhiM(getattr(alljets[j1top], "pt%s"%jesLabel),alljets[j1top].p4().Eta(), alljets[j1top].p4().Phi(), alljets[j1top].p4().M())
            top2 = ROOT.TLorentzVector(); top2.SetPtEtaPhiM(getattr(alljets[j2top], "pt%s"%jesLabel),alljets[j2top].p4().Eta(), alljets[j2top].p4().Phi(), alljets[j2top].p4().M())
            top3 = ROOT.TLorentzVector(); top3.SetPtEtaPhiM(getattr(alljets[j3top], "pt%s"%jesLabel),alljets[j3top].p4().Eta(), alljets[j3top].p4().Phi(), alljets[j3top].p4().M())
            HadTop = top1+top2+top3
        return HadTop

    def analyze(self, event):

        # Some useful input parameters
        year=getattr(event,'year')
        btagvetoval=HiggsRecoTTHbtagwps['DeepFlav_%d_%s'%(year,self.btagDeepCSVveto)][1]

        nAllLeps = event.nLepGood
        nRecleanedLeps = event.nLepFO_Recl
        recleanedLepsIdxs = event.iLepFO_Recl
        allLeps = Collection(event,"LepGood","nLepGood")
        leps = [allLeps[recleanedLepsIdxs[i]] for i in xrange(nRecleanedLeps)]
        alljets = [x for x in Collection(event,"JetSel_Recl","nJetSel_Recl")]

        try:
            #genpar = Collection(event,"GenPartFlav","nGenPartFlav")
            gen = leps[0].genPartFlav
            isData = False
        except:
            isData = True
        
        if(isData): self.systsJEC = {0:""}

        (met, met_phi)  = event.MET_pt, event.MET_phi # what about propagation of JES to MET?

        for jesLabel in self.systsJEC.values():

            if len(leps) < 2                      :
                self.setDefault(event, jesLabel)
                continue

            # Build the jets
            jets = []
            for j in alljets:
                if alljets.index(j) in [int(getattr(event,"BDThttTT_eventReco_iJetSel1%s"%jesLabel)),int(getattr(event,"BDThttTT_eventReco_iJetSel2%s"%jesLabel)),int(getattr(event,"BDThttTT_eventReco_iJetSel3%s"%jesLabel))]:
                    setattr(j, 'fromHadTop', True)
                else:
                    setattr(j, 'fromHadTop', False)
                if(getattr(j, "pt%s"%jesLabel) < 25): continue

                jets.append(j)
            # Store all jets
            jet_pts=[]; jet_etas=[]; jet_phis=[]; jet_masses=[]; jet_btagdiscrs=[]; jet_ishadtops=[]
            def my_sort(j):
                #return j.btagDeepFlavB
                return getattr(j, "pt%s"%jesLabel)
            jets.sort(key=my_sort, reverse=True)
            #if(len(jets) > 5): continue
            for j in jets:
                jet_pts.append(getattr(j, "pt%s"%jesLabel))
                jet_etas.append(j.eta)
                jet_phis.append(j.phi)
                jet_masses.append(j.mass)
                jet_btagdiscrs.append( j.btagDeepFlavB > btagvetoval) 
                jet_ishadtops.append(j.fromHadTop)
            #if(len(jets) < 2): 
            #    self.setDefault(event, jesLabel)
            #    continue


            all5_jets = TLorentzVector()
            all5_jets.SetPtEtaPhiM(0,0,0,0)

            #all_jet_feature.sort()

            for i in range(self.njet):

                if(i < len(jets)): 
                    j5_t = TLorentzVector(0,0,0,0)
                    j5_t.SetPtEtaPhiM(getattr(jets[i], "pt%s"%jesLabel), jets[i].eta, jets[i].phi, jets[i].mass)
                    all5_jets = all5_jets + j5_t
 
            score = getattr(event,"BDThttTT_eventReco_mvaValue%s"%jesLabel)
            
            HadTop = self.buildHadronicTop(event, score, alljets, jesLabel)


            
            #Compute met from events
            all_jets = TLorentzVector()
            all_jets.SetPtEtaPhiM(0,0,0,0)
            all_jets_lep = TLorentzVector()
            all_jets_lep.SetPtEtaPhiM(0,0,0,0)
            for j in range(len(jets)):
                jp = TLorentzVector(0,0,0,0)
                jp.SetPtEtaPhiM(getattr(jets[j], "pt%s"%jesLabel), jets[j].eta, jets[j].phi, jets[j].mass)
                all_jets = all_jets + jp
       
            all_jets_sum = all_jets
            all_leps = TLorentzVector()
            all_leps.SetPtEtaPhiM(0,0,0,0)
            for l in range(len(leps)):
                all_leps = all_leps + leps[l].p4()
                #                all_jets_lep = all_jets_sum + leps[l].p4()

                #            for iLep in range(self.nlep):
            for iLep in range(self.nlep):
                part = leps[iLep].p4()
                self.out.fillBranch('%sLep%s%s_pt'  %(self.label,iLep,jesLabel), part.Pt()  )
                self.out.fillBranch('%sLep%s%s_eta' %(self.label,iLep,jesLabel), part.Eta() )
                self.out.fillBranch('%sLep%s%s_phi' %(self.label,iLep,jesLabel), part.Phi() )
                self.out.fillBranch('%sLep%s%s_mass'%(self.label,iLep,jesLabel), part.M()   )
                
            higgs_mass = TLorentzVector()
            _met = TLorentzVector()
            _met.SetPtEtaPhiM(event.MET_pt, 0, event.MET_phi, 0)
            higgs_reco_mass = all_jets_lep + _met + all_jets
            
            #self.out.fillBranch('%shiggs_reco_mass%s'%(self.label,jesLabel), higgs_reco_mass.M())

            #Compute met from events
            more5_jets = TLorentzVector()
            more5_jets.SetPtEtaPhiM(0,0,0,0)
            if(len(jets) > self.njet):
                #more5_jets = TLorentzVector()
                #more5_jets.SetPtEtaPhiM(0,0,0,0)
                for j in range(self.njet, len(jets)):
                    m5 = TLorentzVector(0,0,0,0)
                    m5.SetPtEtaPhiM(getattr(jets[j], "pt%s"%jesLabel), jets[j].eta, jets[j].phi, jets[j].mass                    )
                    more5_jets = more5_jets + m5


#            self.out.fillBranch('%sHTXS_Higgs_pt'%(self.label), getattr(event,"HTXS_Higgs_pt") if(isData == False) else -99)


            # I must patch these two to fill only for TTH, otherwise the friend does not exist etc. Maybe produce friend also for background
            #self.out.fillBranch('%sHgen_vis_pt%s'  %(self.label,jesLabel), getattr(event,'Hreco_pTTrueGen'))
            #self.out.fillBranch('%sHgen_tru_pt%s'  %(self.label,jesLabel), getattr(event,'Hreco_pTTrueGenPlusNu')) # the same as HTXS_Higgs_pt
            #dnn_pred = self.model_regression.predict(np.transpose(np.array([ [leps[0].p4().Pt()], [leps[0].p4().Eta()], [leps[0].p4().Phi()], [leps[1].p4().Pt()], [leps[1].p4().Eta()], [leps[1].p4().Phi()], [HadTop.Pt()], [HadTop.Eta()], [HadTop.Phi()], [score], [met], [all_jets.Pt()], [all_jets.Eta()], [all_jets.Phi()], [more5_jets.Pt()], [more5_jets.Eta()], [more5_jets.Phi()], [all5_jets.Pt()], [all5_jets.Eta()], [all5_jets.Phi()], [met_phi] ])))
            #print len(leps)
            if(len(leps) >=2 ):
                dnn_pred = self.model_regression.predict(np.transpose(np.array([ [leps[0].p4().Pt()], [leps[0].p4().Eta()], [leps[0].p4().Phi()], [leps[1].p4().Pt()], [leps[1].p4().Eta()], [leps[1].p4().Phi()], [HadTop.Pt() if HadTop is True else -99], [HadTop.Eta() if HadTop is True else -99], [HadTop.Phi() if HadTop is True else -99], [score], [met], [all_jets_lep.Pt()], [all_jets_lep.Eta()], [all_jets_lep.Phi()], [more5_jets.Pt() if(len(jets)) > 5 else -99], [more5_jets.Eta() if(len(jets)) > 5 else -99], [more5_jets.Phi() if(len(jets)) > 5 else -99], [all5_jets.Pt() if(len(jets)) <= 5 else -99], [all5_jets.Eta() if(len(jets)) <= 5 else -99], [all5_jets.Phi() if(len(jets)) <= 5 else -99], [met_phi] ])))
                self.out.fillBranch('%shiggs_reco_mass%s'%(self.label,jesLabel), higgs_reco_mass.M())
                self.out.fillBranch('%sHTXS_Higgs_pt'%(self.label), getattr(event,"HTXS_Higgs_pt") if(isData == False) else -99)
            else:
                dnn_pred = -99.
                self.out.fillBranch('%shiggs_reco_mass%s'%(self.label,jesLabel), -99.)
                self.out.fillBranch('%sHTXS_Higgs_pt'%(self.label), -99.)
            #print dnn_pred
            self.out.fillBranch('%sdnn_prediction%s'%(self.label,jesLabel), dnn_pred)
            
        return True

ttH_2lss_dnn_pt_regression_EFT = lambda : Class_ttH_2lss_dnn_pt_regression_EFT(label='Hreco_',
                                                         btagDeepCSVveto = 'M')
