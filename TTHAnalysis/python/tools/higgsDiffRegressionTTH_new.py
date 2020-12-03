from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from PhysicsTools.Heppy.physicsobjects.Jet import _btagWPs as HiggsRecoTTHbtagwps

import ROOT, itertools
from ROOT import *
import numpy as np

class HiggsDiffRegressionTTH(Module):
    def __init__(self,label="_Recl", cut_BDT_rTT_score = 0.0, btagDeepCSVveto = 'M', doSystJEC=False):
        self.label = label
        self.cut_BDT_rTT_score = cut_BDT_rTT_score
        self.btagDeepCSVveto = btagDeepCSVveto
        self.branches = []
        self.systsJEC = {0:"", 1:"_jesTotalCorrUp", -1:"_jesTotalCorrDown", 2:"_jesTotalUnCorrUp", -2:"_jesTotalUnCorrDown"} if doSystJEC else {0:""}
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Independent on JES

        # Somehow dependent on JES

        for jesLabel in self.systsJEC.values():
            
            # Counters
            self.out.branch('%snLeps%s'%(self.label,jesLabel), 'I') 
            self.out.branch('%snJets%s'%(self.label,jesLabel), 'I')

            # Leptons and the precomputed hadronic top
            for suffix in ["_pt", "_eta", "_phi", "_mass"]:
                for iLep in range(2):
                    self.out.branch('%sLep%s%s%s'%(self.label,iLep,jesLabel,suffix)   , 'F') 
                self.out.branch('%sHadTop%s%s'%(self.label,jesLabel,suffix), 'F')

            # Jets
            for suffix in ["_pt", "_eta", "_phi", "_mass", "_isbtagged", "_ishadtop"]:
                self.out.branch('%sJet%s%s'%(self.label,jesLabel,suffix), 'F' , 24, '%snJets%s'%(self.label,jesLabel))
                for iJet in range(5):
                    self.out.branch('%sJet%s%s%s'%(self.label,iJet,jesLabel,suffix)   , 'F')

            self.out.branch('%sTopScore%s'%(self.label,jesLabel)      , 'F')      
            self.out.branch('%smet%s'%(self.label,jesLabel)           , 'F')       
            self.out.branch('%smet_phi%s'%(self.label,jesLabel)       , 'F')
            self.out.branch('%sHTXS_Higgs%s_pt'%(self.label,jesLabel) , 'F')
            self.out.branch('%sHTXS_Higgs%s_y'%(self.label,jesLabel)  , 'F')
            self.out.branch('%sHgen_vis_pt%s'%(self.label,jesLabel)   , 'F')
            self.out.branch('%sHgen_tru_pt%s'%(self.label,jesLabel)   , 'F')

            self.out.branch('%snNuFromWFromH%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromH_Pt%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromH_Eta%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromH_Phi%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromH_M%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%snNuFromWFromT%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%snNuFromWFromT_lenpt%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromT_Pt%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromT_Eta%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromT_Phi%s'%(self.label,jesLabel)        , 'F')
            self.out.branch('%sNuFromWFromT_M%s'%(self.label,jesLabel)        , 'F')


            self.out.branch('%sMore5_Jets_pt'%(self.label), 'F'   )
            self.out.branch('%sMore5_Jets_eta'%(self.label), 'F'   )
            self.out.branch('%sMore5_Jets_phi'%(self.label), 'F'   )
            self.out.branch('%sMore5_Jets_mass'%(self.label), 'F'   )

            self.out.branch('%sAll5_Jets_pt'%(self.label), 'F'   )
            self.out.branch('%sAll5_Jets_eta'%(self.label), 'F'   )
            self.out.branch('%sAll5_Jets_phi'%(self.label), 'F'   )
            self.out.branch('%sAll5_Jets_mass'%(self.label), 'F'   )


            self.out.branch('%sJets_plus_Lep_pt'%(self.label), 'F'   )
            self.out.branch('%sJets_plus_Lep_eta'%(self.label), 'F'   )
            self.out.branch('%sJets_plus_Lep_phi'%(self.label), 'F'   )
            self.out.branch('%sJets_plus_Lep_mass'%(self.label), 'F'   )
            self.out.branch('%sMet_calc_pt'%(self.label), 'F'   )
            self.out.branch('%sMet_calc_eta'%(self.label), 'F'   )
            self.out.branch('%sMet_calc_phi'%(self.label), 'F'   )
            self.out.branch('%sMet_calc_mass'%(self.label), 'F'   )

            self.out.branch('%sevt_tag%s'%(self.label,jesLabel)       , 'F')       
            
            # Tell the network how to connect the lepton and the closest jet
            for iLep in range(2):
                self.out.branch('%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_pt%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel) , 'F')
                self.out.branch('%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel) , 'F')

            # In principle I need only the DeltaRl0l1 and the previously defined "closest" vars
            for var in ['DeltaRl0l1',
                        #'DeltaRl0j0', 'DeltaRl0j1', 'DeltaRl0j2', 'DeltaRl0j3', 'DeltaRl0j4', 'DeltaRl0j5', 'DeltaRl0j6', 
                        #'DeltaRl1j0', 'DeltaRl1j1', 'DeltaRl1j2', 'DeltaRl1j3', 'DeltaRl1j4', 'DeltaRl1j5', 'DeltaRl1j6', 
                        #'DeltaRj0j0', 'DeltaRj0j1', 'DeltaRj0j2', 'DeltaRj0j3', 'DeltaRj0j4', 'DeltaRj0j5', 'DeltaRj0j6', 
                        #'DeltaRj1j0', 'DeltaRj1j1', 'DeltaRj1j2', 'DeltaRj1j3', 'DeltaRj1j4', 'DeltaRj1j5', 'DeltaRj1j6', 
                        #'DeltaRj2j0', 'DeltaRj2j1', 'DeltaRj2j2', 'DeltaRj2j3', 'DeltaRj2j4', 'DeltaRj2j5', 'DeltaRj2j6', 
                        #'DeltaRj3j0', 'DeltaRj3j1', 'DeltaRj3j2', 'DeltaRj3j3', 'DeltaRj3j4', 'DeltaRj3j5', 'DeltaRj3j6', 
                        #'DeltaRj4j0', 'DeltaRj4j1', 'DeltaRj4j2', 'DeltaRj4j3', 'DeltaRj4j4', 'DeltaRj4j5', 'DeltaRj4j6', 
                        #'DeltaRj5j0', 'DeltaRj5j1', 'DeltaRj5j2', 'DeltaRj5j3', 'DeltaRj5j4', 'DeltaRj5j5', 'DeltaRj5j6', 
                        #'DeltaRj6j0', 'DeltaRj6j1', 'DeltaRj6j2', 'DeltaRj6j3', 'DeltaRj6j4', 'DeltaRj6j5', 'DeltaRj6j6', 
                    ]:
                self.out.branch('%s%s%s'%(self.label,var,jesLabel), 'F')

    def buildHadronicTop(self, event, score, alljets, jesLabel):
        HadTop=None
        if score>self.cut_BDT_rTT_score:
            j1top = int(getattr(event,"BDThttTT_eventReco_iJetSel1%s"%jesLabel))
            j2top = int(getattr(event,"BDThttTT_eventReco_iJetSel2%s"%jesLabel))
            j3top = int(getattr(event,"BDThttTT_eventReco_iJetSel3%s"%jesLabel))
            # Build hadronic top
            top1 = ROOT.TLorentzVector(); top1.SetPtEtaPhiM(alljets[j1top].p4().Pt(),alljets[j1top].p4().Eta(), alljets[j1top].p4().Phi(), alljets[j1top].p4().M())
            top2 = ROOT.TLorentzVector(); top2.SetPtEtaPhiM(alljets[j2top].p4().Pt(),alljets[j2top].p4().Eta(), alljets[j2top].p4().Phi(), alljets[j2top].p4().M())
            top3 = ROOT.TLorentzVector(); top3.SetPtEtaPhiM(alljets[j3top].p4().Pt(),alljets[j3top].p4().Eta(), alljets[j3top].p4().Phi(), alljets[j3top].p4().M())
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

        # Enforce selection (non-jet part)
        if len(leps) < 2                      : return False
        if leps[0].pt < 25 or leps[1].pt < 15 : return False # pt2515
        if event.nLepTight_Recl > 2           : return False # exclusive
        if leps[0].pdgId*leps[1].pdgId < 0    : return False # same-sign
        if abs(event.mZ1_Recl-91.2)<10        : return False # Z_veto
        
        # Make sure we have prompt leptons
        if leps[0].genPartFlav != 1 and leps[0].genPartFlav != 15 : return False
        if leps[1].genPartFlav != 1 and leps[1].genPartFlav != 15 : return False
        
        # No final states with taus, for the moment
        if event.nTauSel_Recl_Tight > 0                           : return False

        (met, met_phi)  = event.MET_pt, event.MET_phi # what about propagation of JES to MET?

        for jesLabel in self.systsJEC.values():

            # Enforce selection (jets part)
            if not ((getattr(event,"nJet25%s_Recl"%jesLabel)>=3 and (getattr(event,"nBJetLoose25%s_Recl"%jesLabel)>= 2 or getattr(event,"nBJetMedium25%s_Recl"%jesLabel)>= 1)) or (getattr(event,"nBJetMedium25%s_Recl"%jesLabel) >= 1 and (getattr(event,"nJet25%s_Recl"%jesLabel)+getattr(event,"nFwdJet%s_Recl"%jesLabel)-getattr(event,"nBJetLoose25%s_Recl"%jesLabel)) > 0)) : continue

            # Build the jets
            jets = []
            for j in alljets:
                if alljets.index(j) in [int(getattr(event,"BDThttTT_eventReco_iJetSel1%s"%jesLabel)),int(getattr(event,"BDThttTT_eventReco_iJetSel2%s"%jesLabel)),int(getattr(event,"BDThttTT_eventReco_iJetSel3%s"%jesLabel))]:
                    setattr(j, 'fromHadTop', True)
                else:
                    setattr(j, 'fromHadTop', False)
                if j.pt < 25: continue
                #if j.btagDeepFlavB < btagvetoval: continue
                jets.append(j)

            # Store all jets
            jet_pts=[]; jet_etas=[]; jet_phis=[]; jet_masses=[]; jet_isbtaggeds=[]; jet_ishadtops=[]
            def my_sort(j):
                return j.btagDeepFlavB
            jets.sort(key=my_sort)
            #if(len(jets) > 5): continue
            for j in jets:
                jet_pts.append(j.pt)
                jet_etas.append(j.eta)
                jet_phis.append(j.phi)
                jet_masses.append(j.mass)
                jet_isbtaggeds.append( j.btagDeepFlavB > btagvetoval) 
                jet_ishadtops.append(j.fromHadTop)

            
            all5_jets = TLorentzVector()
            all5_jets.SetPtEtaPhiM(0,0,0,0)

            #all_jet_feature.sort()
            self.out.fillBranch('%snJets%s'%(self.label,jesLabel)        , len(jets))
            for i in range(5):
                #part = jet_isbtaggeds[i]
                #self.out.fillBranch('%snJets%s'%(self.label,i,jesLabel)        , len(jets))  
                #if i<len(jets) else None

                if(i < len(jets)): all5_jets = all5_jets + jets[i].p4()


                self.out.fillBranch('%sJet%s%s_pt'%(self.label,i,jesLabel)       , jet_pts[i] if i < len(jets)and jet_ishadtops[i] is True else -99)
                self.out.fillBranch('%sJet%s%s_eta'%(self.label,i,jesLabel)      , jet_etas[i] if i < len(jets) and jet_ishadtops[i] is True else -99)
                self.out.fillBranch('%sJet%s%s_phi'%(self.label,i,jesLabel)      , jet_phis[i] if i < len(jets) and jet_ishadtops[i] is True else -99) 
                self.out.fillBranch('%sJet%s%s_mass'%(self.label,i,jesLabel)     , jet_masses[i] if i < len(jets) and jet_ishadtops[i] is True else -99)
                self.out.fillBranch('%sJet%s%s_isbtagged'%(self.label,i,jesLabel), jet_isbtaggeds[i] if i < len(jets) and jet_ishadtops[i] is True else -99)
                self.out.fillBranch('%sJet%s%s_ishadtop'%(self.label,i,jesLabel) , jet_ishadtops[i] if i < len(jets) and jet_ishadtops[i] is True else -99) 

            self.out.fillBranch('%sAll5_Jets_pt'%(self.label), all5_jets.Pt())
            self.out.fillBranch('%sAll5_Jets_eta'%(self.label), all5_jets.Eta())
            self.out.fillBranch('%sAll5_Jets_phi'%(self.label), all5_jets.Phi())
            self.out.fillBranch('%sAll5_Jets_mass'%(self.label), all5_jets.M())


            score = getattr(event,"BDThttTT_eventReco_mvaValue%s"%jesLabel)
            
            HadTop = self.buildHadronicTop(event, score, alljets, jesLabel)

            self.out.fillBranch('%sHadTop%s_pt'  %(self.label,jesLabel), HadTop.Pt()  if HadTop else -99.)
            self.out.fillBranch('%sHadTop%s_eta' %(self.label,jesLabel), HadTop.Eta() if HadTop else -99.)
            self.out.fillBranch('%sHadTop%s_phi' %(self.label,jesLabel), HadTop.Phi() if HadTop else -99.)
            self.out.fillBranch('%sHadTop%s_mass'%(self.label,jesLabel), HadTop.M()   if HadTop else -99.)
            self.out.fillBranch('%sTopScore%s'   %(self.label,jesLabel), score                           ) # by not filling it with -99, a network should be able to learn self.cut_BDT_rTT_score

            # Lepton observables
            self.out.fillBranch('%snLeps%s' %(self.label,jesLabel), len(leps))
            self.out.fillBranch('%sevt_tag%s'%(self.label,jesLabel), leps[0].pdgId*leps[1].pdgId)
            self.out.fillBranch('%sDeltaRl0l1%s' %(self.label,jesLabel), leps[0].p4().DeltaR(leps[1].p4()) if len(leps)>=2 else -99.)

            for iLep in range(2):
                part = leps[iLep].p4()
                self.out.fillBranch('%sLep%s%s_pt'  %(self.label,iLep,jesLabel), part.Pt()  )
                self.out.fillBranch('%sLep%s%s_eta' %(self.label,iLep,jesLabel), part.Eta() )
                self.out.fillBranch('%sLep%s%s_phi' %(self.label,iLep,jesLabel), part.Phi() )
                self.out.fillBranch('%sLep%s%s_mass'%(self.label,iLep,jesLabel), part.M()   )

            #Compute met from events
            all_jets = TLorentzVector()
            all_jets.SetPtEtaPhiM(0,0,0,0)
            for j in range(len(jets)):
                all_jets = all_jets + jets[j].p4()
            
            for l in range(len(leps)):
                all_jets=all_jets + leps[l].p4()

            self.out.fillBranch('%sJets_plus_Lep_pt'%(self.label), all_jets.Pt()   )
            self.out.fillBranch('%sJets_plus_Lep_eta'%(self.label), all_jets.Eta()   )
            self.out.fillBranch('%sJets_plus_Lep_phi'%(self.label), all_jets.Phi()   )
            self.out.fillBranch('%sJets_plus_Lep_mass'%(self.label), all_jets.M()   )
            

            #Compute met from events
            if(len(jets) > 5):
                more5_jets = TLorentzVector()
                more5_jets.SetPtEtaPhiM(0,0,0,0)
                for j in range(5, len(jets)):
                    more5_jets = more5_jets + jets[j].p4()

                self.out.fillBranch('%sMore5_Jets_pt'%(self.label), more5_jets.Pt() if(len(jets)) >=5 else -10   )
                self.out.fillBranch('%sMore5_Jets_eta'%(self.label), more5_jets.Eta() if(len(jets)) >=5 else -10  )
                self.out.fillBranch('%sMore5_Jets_phi'%(self.label), more5_jets.Phi() if(len(jets)) >=5 else -10  )
                self.out.fillBranch('%sMore5_Jets_mass'%(self.label), more5_jets.M() if(len(jets)) >=5 else -10  )
            

            met_calc = TLorentzVector()
            met_calc.SetPtEtaPhiM(met,-all_jets.Eta(),met_phi,0)

            self.out.fillBranch('%sMet_calc_pt'%(self.label), met_calc.Pt()   )
            self.out.fillBranch('%sMet_calc_eta'%(self.label), met_calc.Eta()   )
            self.out.fillBranch('%sMet_calc_phi'%(self.label), met_calc.Phi()   )
            self.out.fillBranch('%sMet_calc_mass'%(self.label), met_calc.M()   )

            # Compute the deltaR and pt of the closest jet to each lepton
            deltaR_closestJet = [999.,999.]
            deltaPt_closestJet = [999.,999.]
            lep_closestJet = []

            for iLep in range(2):
                lp4=leps[iLep].p4()
                for j, jp4 in [(ix,x.p4()) for ix,x in enumerate(jets)]:
                    dr_lj = lp4.DeltaR(jp4)
                    dpt_lj = lp4.Pt()-jp4.Pt()

                    if dr_lj < deltaR_closestJet:
                        deltaR_closestJet[iLep] = dr_lj
                        deltaPt_closestJet[iLep] = dpt_lj
                        lep_closestJet.append(lp4+jp4)

            #print(len(lep_closestJet))
            #if(len(lep_closestJet) == 0): print 'cacca ', deltaR_closestJet[0], deltaR_closestJet[1]
            for iLep in range(2): 
                self.out.fillBranch('%sDeltaRClosestJetToLep%s%s'%(self.label,iLep,jesLabel) ,  deltaR_closestJet[iLep])
                self.out.fillBranch('%sDeltaPtClosestJetToLep%s%s'%(self.label,iLep,jesLabel) , deltaPt_closestJet[iLep])
                self.out.fillBranch('%sClosestJetToLep_pt%s%s'%(self.label,iLep,jesLabel) , lep_closestJet[iLep].Pt() if iLep < len(lep_closestJet) else -99)            # Fill MET and GEN observables        
                self.out.fillBranch('%sClosestJetToLep_eta%s%s'%(self.label,iLep,jesLabel) , lep_closestJet[iLep].Eta() if iLep < len(lep_closestJet) else -99)
                self.out.fillBranch('%sClosestJetToLep_phi%s%s'%(self.label,iLep,jesLabel) , lep_closestJet[iLep].Phi() if iLep < len(lep_closestJet) else -99)
                self.out.fillBranch('%sClosestJetToLep_mass%s%s'%(self.label,iLep,jesLabel) , lep_closestJet[iLep].M() if iLep < len(lep_closestJet) else -99)

            self.out.fillBranch('%smet%s'     %(self.label,jesLabel), met                                ) 
            self.out.fillBranch('%smet_phi%s' %(self.label,jesLabel), met_phi                            )
            self.out.fillBranch('%sHTXS_Higgs_pt%s'%(self.label,jesLabel), getattr(event,"HTXS_Higgs_pt"))
            self.out.fillBranch('%sHTXS_Higgs_y%s' %(self.label,jesLabel), getattr(event,"HTXS_Higgs_y") )
            # I must patch these two to fill only for TTH, otherwise the friend does not exist etc. Maybe produce friend also for background
            self.out.fillBranch('%sHgen_vis_pt%s'  %(self.label,jesLabel), getattr(event,'Hreco_pTTrueGen'))
            self.out.fillBranch('%sHgen_tru_pt%s'  %(self.label,jesLabel), getattr(event,'Hreco_pTTrueGenPlusNu')) # the same as HTXS_Higgs_pt

            self.out.fillBranch('%snNuFromWFromH%s'   %(self.label,jesLabel), getattr(event, 'Hreco_nNuFromWFromH'))
            self.out.fillBranch('%sNuFromWFromH_Pt%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromH_pt')[0])
            self.out.fillBranch('%sNuFromWFromH_Eta%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromH_eta')[0])
            self.out.fillBranch('%sNuFromWFromH_Phi%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromH_phi')[0])
            self.out.fillBranch('%sNuFromWFromH_M%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromH_mass')[0])


            self.out.fillBranch('%snNuFromWFromT_lenpt%s'   %(self.label,jesLabel), len(getattr(event, 'Hreco_NuFromWFromT_pt')))
            self.out.fillBranch('%snNuFromWFromT%s'   %(self.label,jesLabel), getattr(event, 'Hreco_nNuFromWFromT'))
            self.out.fillBranch('%sNuFromWFromT_Pt%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromT_pt')[0])
            self.out.fillBranch('%sNuFromWFromT_Eta%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromT_eta')[0])
            self.out.fillBranch('%sNuFromWFromT_Phi%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromT_phi')[0])
            self.out.fillBranch('%sNuFromWFromT_M%s'   %(self.label,jesLabel), getattr(event, 'Hreco_NuFromWFromT_mass')[0])

        return True

higgsDiffRegressionTTH_new = lambda : HiggsDiffRegressionTTH(label='Hreco_',
                                                         btagDeepCSVveto = 'M')
