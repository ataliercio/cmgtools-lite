### First prepare plots

# From: python ttH-multilepton/ttH_plots.py 2lss_systcheck 2lss_SR --plotmode norm 
# with tweaks for mca and vars
#python mcPlots.py -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skimOnlyMC_v6 --pdir 2lss_systcheck/2lss_SR --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6  --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/2lss_tight.txt  -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],2)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],2)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],2)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^lep(3|4)_.*' --xP '^(3|4)lep_.*' --xP 'kinMVA_3l_.*'    --plotmode norm --sP kinMVA_2lss_ttV_withHj --sP kinMVA_2lss_ttbar_withBDTv8 -j 8
#python mcPlots.py -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skimOnlyMC_v6 --pdir 2lss_systcheck_bins1/2lss_SR --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6  --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/2lss_tight.txt  -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],2)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],2)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],2)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^lep(3|4)_.*' --xP '^(3|4)lep_.*' --xP 'kinMVA_3l_.*'    --plotmode norm --sP kinMVA_bins1_2lss_ttV_withHj --sP kinMVA_bins1_2lss_ttbar_withBDTv8 -j 8


# From: python ttH-multilepton/ttH_plots.py 2lss_systcheck 3l_SR --plotmode norm 
# with tweaks for mca and vars
#python mcPlots.py --pdir 2lss_systcheck/3l_SR --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6 -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skimOnlyMC_v6 --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/3l_tight.txt  -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[2]],LepGood_pt[iLepFO_Recl[2]],LepGood_eta[iLepFO_Recl[2]],3)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],3)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^(2|4)lep_.*' --xP '^lep4_.*' --xP 'kinMVA_2lss_.*'    --plotmode norm  --sP kinMVA_3l_ttbar --sP kinMVA_3l_ttV -j 8 
#python mcPlots.py --pdir 2lss_systcheck_bins1/3l_SR --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6 -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skimOnlyMC_v6 --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/3l_tight.txt  -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[2]],LepGood_pt[iLepFO_Recl[2]],LepGood_eta[iLepFO_Recl[2]],3)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],3)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^(2|4)lep_.*' --xP '^lep4_.*' --xP 'kinMVA_2lss_.*'    --plotmode norm  --sP kinMVA_bins1_3l_ttbar --sP kinMVA_bins1_3l_ttV -j 8 

# From: python ttH-multilepton/ttH_plots.py 2lss_systcheck 3l_SR_prescale --plotmode norm 
# with tweaks for mca and vars
#python mcPlots.py --pdir 2lss_systcheck/3l_SR_prescale --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6 -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skim3l2j2b1B_v6 --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_withMEM_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/3l_tight.txt -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[2]],LepGood_pt[iLepFO_Recl[2]],LepGood_eta[iLepFO_Recl[2]],3)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],3)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^(2|4)lep_.*' --xP '^lep4_.*' --xP 'kinMVA_2lss_.*'  --Fs {P}/8_MEM_v6 --sP kinMVA_3l_ttbar --sP kinMVA_3l_ttV_withMEM --fitRatio 1 --plotmode norm 
#python mcPlots.py --pdir 2lss_systcheck_bins1/3l_SR_prescale --Fs {P}/1_recleaner_230217_v6 --Fs {P}/5_triggerDecision_230217_v6 --Fs {P}/6_bTagSF_v6 --Fs {P}/7_tauTightSel_v6 -P /afs/cern.ch/work/p/peruzzi/tthtrees/TREES_TTH_250117_Summer16_JECV3_noClean_qgV2_skim3l2j2b1B_v6 --Fs {P}/4_BDTv8_Hj_230217_v6 --Fs {P}/3_kinMVA_BDTv8_withMEM_230217_v6  --Fs {P}/2_eventVars_230217_v6 -f -j 8 -l 35.9 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/ttH_2lss3l_triggerdefs.txt --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0 2 --fixRatioRange  --showMCError --rebin 4 --xP 'nT_.*' --xP 'debug_.*' ttH-multilepton/mca-lhe-vars.txt ttH-multilepton/3l_tight.txt -W 'puw2016_nTrueInt_36fb(nTrueInt)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],3)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[2]],LepGood_pt[iLepFO_Recl[2]],LepGood_eta[iLepFO_Recl[2]],3)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],3)*eventBTagSF' ttH-multilepton/2lss_3l_plots.txt --xP '^(2|4)lep_.*' --xP '^lep4_.*' --xP 'kinMVA_2lss_.*'  --Fs {P}/8_MEM_v6 --sP kinMVA_bins1_3l_ttbar --sP kinMVA_bins1_3l_ttV_withMEM --fitRatio 1 --plotmode norm 


# Normalization + shape
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck\",\"2lss_SR\",true,true,true,false\)
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck\",\"3l_SR_prescale\",true,true,true,false\)
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck\",\"3l_SR\",true,true,true,false\)

# Shape only
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck_normscale\",\"2lss_SR\",true,true,true,true\)
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck_normscale\",\"3l_SR_prescale\",true,true,true,true\)
root -l -b ttH-multilepton/makeBdtScalePDFVariations.C\(\"2lss_systcheck_normscale\",\"3l_SR\",true,true,true,true\)


