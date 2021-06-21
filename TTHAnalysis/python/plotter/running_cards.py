import os 
import argparse

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('YES', 'True', 'yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('NO', 'False', 'no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description="prog [options] mc.txt cuts.txt var bins")
parser.add_argument("--year",dest="year",type=str, default=None, help="year")
parser.add_argument("--in_dir",dest="in_dir",type=str, default=None, help="name of input directory, with new variables stored")
parser.add_argument("--out_dir",dest="out_dir",type=str, default=None, help="output directory")
parser.add_argument("--reco_range",dest="range_",type=str, default=None, help="range reco")
parser.add_argument("--txt_cut",dest="cuts",type=str, default=None, help="name of the cut file")
parser.add_argument("--unf_obs",dest="obs",type=str, default=None, help="name of the branch of the observable to unfold")
parser.add_argument("--final_state",dest="fs",type=str, default=None, help="final state, choose one: 2lss, 3l, 2lss1tau or all")
parser.add_argument("--cluster",dest="cluster",type=str2bool, default=None, help="run on cluster")

args = parser.parse_args()

year_ = args.year
reco_range_ = args.range_# "[0,30,60,90,120,160,200,250,300,400,450,550,750,1000]"
in_dir_ = args.in_dir
output = args.out_dir
cuts_ = args.cuts
obs_ = args.obs
fs_ = args.fs
cluster_ = args.cluster

if('2lss' in fs_):
    _2lss = True
    _3l = False
    _2lss1tau = False
if('3l' in fs_):
    _3l = True
    _2lss = False
    _2lss1tau = False
if('2lss1tau' in fs_):
    _2lss1tau = True
    _2lss = False
    _3l = False
if('all' in fs_):
    _2lss1tau = True
    _2lss = True
    _3l = True

#controls on inputs
if('2lss' not in fs_ and '3l' not in fs_ and '2lss1tau' not in fs_ and 'all' not in fs_):
    raise RuntimeError("Error in 'final state', please select a valid argument: 2lss, 3l, 2lss1tau, all")

if('2lss' in str(cuts_) or '3l' in str(cuts_) or '2lss1tau' in str(cuts_)):
    raise RuntimeError("Error with 'cuts', please insert only the postfix of the txt --es 2lss_tight.txt -> --txt_cut tight.txt")
if('.txt' not in str(cuts_)):
    raise RuntimeError("Error with 'cuts', please insert '.txt' at the end of the string")
print year_

if('2016' not in str(year_) and '2017' not in str(year_) and '2018' not in str(year_)):
    raise RuntimeError("Error with 'year', please insert a valid argument: 2016 or 2017 or 2018")


if('2016' in str(year_)):
    luminosity = 35.9
    input_dir = in_dir_+year_
    #input_dir_3l = "/nfs/user/atalier/skim_3l_cp_new/"
    #blah = "{P}/5_BDThtt_reco_new_blah"
if('2017' in str(year_)):
    luminosity = 41.5
    input_dir = in_dir_+year_
    #input_dir_3l = "/nfs/user/atalier/skim_3l_2017_new/"
    #blah = "{P}/5_BDThtt_reco_new_blah"
if('2018' in str(year_)):
    luminosity = 59.7
    input_dir = in_dir_+year_
    #input_dir_3l = "/nfs/user/atalier/skim_3l_2018/"
    #blah = "/nfs/user/atalier/big_ntuples/2018/5_BDThtt_reco_new_blah"


#if('2016' in str(year[0])):
if(_2lss is True):
    print('2lss cards...')
    if(cluster_ is True):
        print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
        os.system('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
    else:
        print('python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6')
        os.system('python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)" --binname ttH_2lss_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6')




if(_3l is True):
    print('3l cards...')
    if(cluster_ is True):
        print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/3l_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
        os.system('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/3l_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
    else:
        print('python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/3l_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6')
        os.system('python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-3l-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/3l_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+'/ --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l  --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L ttH-multilepton/functionsTTH.cc --mcc ttH-multilepton/lepchoice-ttH-FO.txt --mcc ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_3l*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 3, year)" --binname ttH_3l_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6')



if(_2lss1tau is True):
    print('2lss1tau cards...')
    if(cluster_ is True):
        print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss1tau_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc  ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+' --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l --Fs {P}/6_mva_2lss1tau --FMCs {P}/6_mva_tauSFs --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L  ttH-multilepton/functionsTTH.cc --mcc  ttH-multilepton/lepchoice-ttH-FO.txt --mcc  ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter  ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)*TauSel_2lss1tau_SF" --binname ttH_2lss1tau_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
        os.system('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss1tau_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc  ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+' --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l --Fs {P}/6_mva_2lss1tau --FMCs {P}/6_mva_tauSFs --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L  ttH-multilepton/functionsTTH.cc --mcc  ttH-multilepton/lepchoice-ttH-FO.txt --mcc  ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter  ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)*TauSel_2lss1tau_SF" --binname ttH_2lss1tau_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
    else:
        print('sbatch -c 16 -p cp3 --wrap \'python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss1tau_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc  ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+' --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l --Fs {P}/6_mva_2lss1tau --FMCs {P}/6_mva_tauSFs --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L  ttH-multilepton/functionsTTH.cc --mcc  ttH-multilepton/lepchoice-ttH-FO.txt --mcc  ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter  ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)*TauSel_2lss1tau_SF" --binname ttH_2lss1tau_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6\'')
        os.system('python makeShapeCardsNew_unfolding.py ttH-multilepton/mca-2lss-mcdata-frdata-diff_signal_fr.txt ttH-multilepton/2lss1tau_'+str(cuts_)+' "'+str(obs_)+'" "'+str(reco_range_)+'" --unc  ttH-multilepton/systsUnc.txt --amc --xu CMS_ttHl_TTZ_lnU,CMS_ttHl_TTW_lnU -P /nfs/user/pvischia/tth/v6/NanoTrees_TTH_090120_091019_v6_skim2lss_forCP/'+str(year_)+' --Fs '+str(input_dir)+' --FMCs {P}/0_jmeUnc_v1 --FDs {P}/1_recl --FMCs {P}/1_recl_allvars --FMCs {P}/2_btag_SFs --FMCs {P}/2_scalefactors_lep_fixed --Fs {P}/3_tauCount --Fs {P}/4_evtVars --FMCs {P}/5_BDThtt_reco_new_blah --Fs {P}/6_mva2lss --Fs {P}/6_mva3l --Fs {P}/6_mva4l --Fs {P}/6_mva_2lss1tau --FMCs {P}/6_mva_tauSFs --xf TTTW --xf TTWH --tree NanoAOD --s2v -j 16 -l '+str(luminosity)+' -f --WA prescaleFromSkim --split-factor=-1  --od '+str(output)+'  -L  ttH-multilepton/functionsTTH.cc --mcc  ttH-multilepton/lepchoice-ttH-FO.txt --mcc  ttH-multilepton/mcc-METFixEE2017.txt --plotgroup data_fakes+=.*_promptsub --neg   --threshold 0.01 --asimov signal --filter  ttH-multilepton/filter-processes.txt  -W "L1PreFiringWeight_Nom*puWeight*btagSF_shape*leptonSF_2lss*triggerSF_ttH(LepGood1_pdgId, LepGood1_conePt, LepGood2_pdgId, LepGood2_conePt, 2, year)*TauSel_2lss1tau_SF" --binname ttH_2lss1tau_'+str(year_)+' --year '+str(year_)+' --tikhonov_unfolding pt_0_60,pt_60_120,pt_120_200,pt_200_300,pt_300_450,pt_450_inf rBin1,rBin2,rBin3,rBin4,rBin5,rBin6')
