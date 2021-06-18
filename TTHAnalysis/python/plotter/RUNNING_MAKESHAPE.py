import os 
import argparse
parser = argparse.ArgumentParser(description="prog [options] mc.txt cuts.txt var bins")
parser.add_argument("--year",dest="year",type=str, default=None, help="year")
#parser.add_argument("--channel",dest="channel",type=str, default=None, help="select one channel")
#parser.add_argument("--all_channel",dest="all_channel",type=str, default=None, help="running on all the channels")
parser.add_argument("--out_dir",dest="out_dir",type=str, default=None, help="output directory")
parser.add_argument("--reco_range",dest="range_",type=str, default=None, help="range")

#(options, args) = parser.parser_args()

args = parser.parse_args()
#if options.list_card:
#card_list = args.list_card.split(",")import argparse


year_ = args.year
reco_range_ = args.range_# "[0,30,60,90,120,160,200,250,300,400,450,550,750,1000]"

_2lss = True
_3l = True
output = args.out_dir#"cards/all_years"


if('2016' in str(year_)):
    luminosity = 35.9
    input_dir_2lss = "/nfs/user/atalier/skim_2lss_cp_new/"
    input_dir_3l = "/nfs/user/atalier/skim_3l_cp_new/"
    blah = "{P}/5_BDThtt_reco_new_blah"
if('2017' in str(year_)):
    luminosity = 41.5
    input_dir_2lss = "/nfs/user/atalier/skim_2lss_2017/"
    input_dir_3l = "/nfs/user/atalier/skim_3l_2017_new/"
    blah = "{P}/5_BDThtt_reco_new_blah"
if('2018' in str(year_)):
    luminosity = 59.7
    input_dir_2lss = "/nfs/user/atalier/skim_2lss_2018/"
    input_dir_3l = "/nfs/user/atalier/skim_3l_2018/"
    blah = "/nfs/user/atalier/big_ntuples/2018/5_BDThtt_reco_new_blah"


#if('2016' in str(year[0])):
if(_2lss is True):
    print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff.txt ttH-multilepton/2lss_tight_legacy.txt "Hreco_dnn_prediction" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir_2lss)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs '+str(blah)+' --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6 \'')

#    os.system('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff.txt ttH-multilepton/2lss_tight_legacy.txt "Hreco_dnn_prediction" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir_2lss)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs '+str(blah)+' --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6 \'')

print '3L' 

if(_3l is True):
    print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff.txt ttH-multilepton/3l_tight_legacy.txt "Hreco_dnn_prediction" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir_3l)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs '+str(blah)+' --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')

#    os.system('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff.txt ttH-multilepton/3l_tight_legacy.txt "Hreco_dnn_prediction" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir_3l)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs '+str(blah)+' --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')

