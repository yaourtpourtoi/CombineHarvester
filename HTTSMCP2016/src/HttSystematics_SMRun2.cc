#include "CombineHarvester/HTTSMCP2016/interface/HttSystematics_SMRun2.h"
#include <vector>
#include <string>
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"

using namespace std;

namespace ch {
    
    using ch::syst::SystMap;
    using ch::syst::SystMapAsymm;
    using ch::syst::era;
    using ch::syst::channel;
    using ch::syst::bin_id;
    using ch::syst::process;
    using ch::syst::bin;
    using ch::JoinStr;
    
    void __attribute__((optimize("O0"))) AddSMRun2Systematics(CombineHarvester & cb, int control_region, bool ttbar_fit, bool no_jec_split) {
        // Create a CombineHarvester clone that only contains the signal
        // categories
        //
        // cb_sig is unused at the moment, (was it ever used in this analysis?) delete?
        //CombineHarvester cb_sig = cb.cp();
        //
        //
        
        
        std::vector<std::string> sig_procs = {"ggH_htt","qqH_htt","WH_htt","ZH_htt","reweighted_ggH_htt_0PM", "reweighted_ggH_htt_0M", "reweighted_ggH_htt_0Mf05ph0","qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt","vbf125_powheg","qqH_sm_htt125", "qqH_ps_htt", "qqH_mm_htt","wh125_powheg","zh125_powheg","WH_sm_htt125","ZH_sm_htt125", "WH_ps_htt","ZH_ps_htt","WH_mm_htt","ZH_mm_htt","WH_sm_htt","ZH_sm_htt","WH_ps_htt","ZH_ps_htt","WH_mm_htt","ZH_mm_htt", "ggHsm_jhu_htt","ggHps_jhu_htt","ggHmm_jhu_htt","ggH_ph_htt" ,"qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt", "ggH_sm_htt", "ggH_ps_htt", "ggH_mm_htt"};
        std::vector<std::string> ggH_sig_procs = {"ggH_htt","reweighted_ggH_htt_0PM", "reweighted_ggH_htt_0M", "reweighted_ggH_htt_0Mf05ph0","ggHsm_jhu_htt","ggHps_jhu_htt","ggHmm_jhu_htt","ggH_ph_htt","ggH_sm_htt", "ggH_ps_htt", "ggH_mm_htt"};
        std::vector<std::string> qqH_sig_procs = {"qqH_htt","qqH_htt125","qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt", "vbf125_powheg","qqH_sm_htt125", "qqH_ps_htt125", "qqH_mm_htt125"};
        
        // N.B. when adding this list of backgrounds to a nuisance, only
        // the backgrounds that are included in the background process
        // defined in MorphingSM2016.cpp are included in the actual DCs
        // This is a list of all MC based backgrounds
        // QCD is explicitly excluded
        std::vector<std::string> all_mc_bkgs = {
            "ZL","ZLL","ZJ","ZTT","TTJ","TTT","TT",
            "W","VV","VVT","VVJ",
            "ggH_hww125","qqH_hww125","EWKZ","Wfakes"};
        std::vector<std::string> all_mc_bkgs_no_W = {
            "ZL","ZLL","ZJ","ZTT","TTJ","TTT","TT",
            "VV","VVT","VVJ",
            "ggH_hww125","qqH_hww125","EWKZ"};
        std::vector<std::string> all_mc_bkgs_no_TTJ = {
            "ZL","ZLL","ZJ","ZTT","TTT","TT",
            "VV","VVT","VVJ",
            "ggH_hww125","qqH_hww125","EWKZ"};
        std::vector<std::string> embed = {"EmbedZTT"};
        std::vector<std::string> real_tau_mc_bkgs = {"ZTT","TTT","TT","VV","VVT"};
            
        //##############################################################################
        //  lumi
        //##############################################################################
        
        // total lumi uncertainty is 2.5% for 2016, 2.3% for 2017, 2.5% for 2018
        //
        // for correlations using https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb


        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV_2016", "lnN", SystMap<>::init(1.022));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_2017", "lnN", SystMap<>::init(1.020));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_13TeV_2018", "lnN", SystMap<>::init(1.015));


        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV_XY", "lnN", SystMap<>::init(1.009));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_XY", "lnN", SystMap<>::init(1.008));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_13TeV_XY", "lnN", SystMap<>::init(1.020));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_LS", "lnN", SystMap<>::init(1.003));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_13TeV_LS", "lnN", SystMap<>::init(1.002));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV_BBD", "lnN", SystMap<>::init(1.004));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_BBD", "lnN", SystMap<>::init(1.004));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV_DB", "lnN", SystMap<>::init(1.005));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_DB", "lnN", SystMap<>::init(1.005));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_BCC", "lnN", SystMap<>::init(1.003));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_13TeV_BCC", "lnN", SystMap<>::init(1.002));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV_GS", "lnN", SystMap<>::init(1.004));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV_GS", "lnN", SystMap<>::init(1.001));


        //##############################################################################
        //  trigger   
        //##############################################################################
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_mt_13TeV", "lnN", SystMap<>::init(1.02));
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_et_13TeV", "lnN", SystMap<>::init(1.02)); 
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs,embed})).channel({"em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_em_13TeV", "lnN", SystMap<>::init(1.04)); 



//        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
//                                             "CMS_eff_trigger_tt_13TeV", "lnN", SystMap<>::init(1.1));

        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM11_13TeV", "shape", SystMap<>::init(1.00));


        // additional uncertainties due to tau SF on cross-triggers - decorrelated for embedded samples in Morphing code
        //  Add back when this is working properly (i.e when you have fixed issue with nan values in weights)
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM11_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM11_13TeV", "shape", SystMap<>::init(1.00));        

        //##############################################################################
        //  Electron, muon and tau Id  efficiencies
        //##############################################################################
        
        cb.cp().AddSyst(cb, "CMS_eff_m", "lnN", SystMap<channel, process>::init
                        ({"mt","mt_2016","mt_2017","mt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, JoinStr({sig_procs, all_mc_bkgs,embed}),  1.02));
        
        // embedded selection efficiency
        cb.cp().AddSyst(cb, "CMS_eff_m_embedsel", "lnN", SystMap<channel, process>::init
                        ({"em","em_2016","em_2017","em_2018","et","et_2016","et_2017","et_2018","tt","tt_2016","tt_2017","tt_2018","mt","mt_2016","mt_2017","mt_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, embed,  1.04)); 
        
        cb.cp().AddSyst(cb, "CMS_eff_e", "lnN", SystMap<channel, process>::init
                        ({"et","et_2016","et_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, JoinStr({sig_procs, all_mc_bkgs,embed}),       1.02));



        // due to different treatments of embedding and MC uncertainties for the tau ID they are included seperatly for now


        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM11_13TeV", "shape", SystMap<>::init(1.00));

        //cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin1_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin3_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin4_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs,embed})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin5_13TeV", "shape", SystMap<>::init(1.00));

        // the uncorrelated part by channel is due to anti-electron and anti-muon disctiminant, it is the same for MC and embedding so decorrelated in the Morphing       
 
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_mt_13TeV", "lnN", SystMap<>::init(1.03));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                             "CMS_eff_t_et_13TeV", "lnN", SystMap<>::init(1.03));

        // 2 real taus
        cb.cp().process(JoinStr({sig_procs, {"ZTT","VVT","TTT","EWKZ"}, embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.06));

        // 1 real tau
        cb.cp().process(JoinStr({{"TTJ","ZJ","VVJ","W","Wfakes","TTfakes"}})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.03));


        
        //##############################################################################
        //  b tag and mistag rate  efficiencies 
        //##############################################################################

        cb.cp().process({"TTT","TT","VVT","VV"}).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).AddSyst(cb,
                                             "CMS_eff_b_13TeV", "shape", SystMap<>::init(1.00));
        
        //##############################################################################
        //  Electron, muon and tau energy Scale
        //##############################################################################
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed/*, {"jetFakes", "QCD"}*/})).channel({"et","et_2016","et_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                             "CMS_scale_e_13TeV", "shape", SystMap<>::init(1.00));
       

 
        // Decay Mode based TES Settings
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_3prong1pizero_13TeV", "shape", SystMap<>::init(1.00));

        // Muon 
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed, /*{"jetFakes", "QCD"}*/})).channel({"mt","mt_2016","mt_2017","mt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                             "CMS_scale_mu_13TeV", "shape", SystMap<>::init(1.00));

        //##############################################################################
        //  Embedded uncertainty on ttbar contamination (and VV contamination)
        //##############################################################################        
        cb.cp().process(embed).AddSyst(cb,"CMS_ttbar_embeded_13TeV", "shape", SystMap<>::init(1.00));
 
        //##############################################################################
        //  jet and met energy Scale
        //##############################################################################


        // MET Systematic shapes - recoil uncertainties for recoil corrected met, unclustered energy uncertainty for samples with no recoil correction, jes uncertainties propogated to met for samples with no recoil correction
        cb.cp().process({"TT","TTJ","TTT","VV","VVJ","VVT"}).channel({"et","et_2016","mt","mt_2016","tt","tt_2016"},false).AddSyst(cb,
                                                  "CMS_scale_met_unclustered_13TeV", "shape", SystMap<>::init(1.00));


        cb.cp().process(JoinStr({sig_procs, {"ZTT","ZLL","ZL","ZJ","EWKZ","W"}})).AddSyst(cb,
                                                  "CMS_htt_boson_reso_met_13TeV", "shape", SystMap<>::init(1.00)); 
        cb.cp().process(JoinStr({sig_procs, {"ZTT","ZLL","ZL","ZJ","EWKZ","W"}})).AddSyst(cb,
                                                  "CMS_htt_boson_scale_met_13TeV", "shape", SystMap<>::init(1.00));
 

      
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_Absolute_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_BBEC1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_EC2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_FlavorQCD_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_HF_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_RelativeBal_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","em","em_2016","et","et_2016","mt","mt_2016"}).AddSyst(cb,"CMS_scale_j_Absolute_2016_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","em_2017","et_2017","mt_2017"}).AddSyst(cb,"CMS_scale_j_Absolute_2017_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","em_2018","et_2018","mt_2018"}).AddSyst(cb,"CMS_scale_j_Absolute_2018_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","em","em_2016","et","et_2016","mt","mt_2016"}).AddSyst(cb,"CMS_scale_j_BBEC1_2016_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","em_2017","et_2017","mt_2017"}).AddSyst(cb,"CMS_scale_j_BBEC1_2017_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","em_2018","et_2018","mt_2018"}).AddSyst(cb,"CMS_scale_j_BBEC1_2018_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","em","em_2016","et","et_2016","mt","mt_2016"}).AddSyst(cb,"CMS_scale_j_EC2_2016_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","em_2017","et_2017","mt_2017"}).AddSyst(cb,"CMS_scale_j_EC2_2017_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","em_2018","et_2018","mt_2018"}).AddSyst(cb,"CMS_scale_j_EC2_2018_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","em","em_2016","et","et_2016","mt","mt_2016"}).AddSyst(cb,"CMS_scale_j_HF_2016_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","em_2017","et_2017","mt_2017"}).AddSyst(cb,"CMS_scale_j_HF_2017_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","em_2018","et_2018","mt_2018"}).AddSyst(cb,"CMS_scale_j_HF_2018_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","em","em_2016","et","et_2016","mt","mt_2016"}).AddSyst(cb,"CMS_scale_j_RelativeSample_2016_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","em_2017","et_2017","mt_2017"}).AddSyst(cb,"CMS_scale_j_RelativeSample_2017_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","em_2018","et_2018","mt_2018"}).AddSyst(cb,"CMS_scale_j_RelativeSample_2018_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_res_j_13TeV", "shape", SystMap<>::init(1.00));
 
        
        //##############################################################################
        //  Background normalization uncertainties
        //##############################################################################
        
        //   Diboson  Normalisation - fully correlated
        cb.cp().process({"VV","VVT","VVJ"}).AddSyst(cb,
                                        "CMS_htt_vvXsec_13TeV", "lnN", SystMap<>::init(1.05));

        cb.cp().process({"ZTT","ZJ","ZL","ZLL"}).AddSyst(cb,
                                        "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.02));        
 
        cb.cp().process({"EWKZ"}).AddSyst(cb,
                                        "CMS_htt_ewkzXsec_13TeV", "lnN", SystMap<>::init(1.05));

        //if (! ttbar_fit){
        //   ttbar Normalisation - fully correlated
	    cb.cp().process({"TT","TTT","TTJ"}).AddSyst(cb,
					  "CMS_htt_tjXsec_13TeV", "lnN", SystMap<>::init(1.042));
        //}

        // W norm, just for em, tt and the mm region where MC norm is from MC
        
        cb.cp().process({"W"}).channel({"em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                                       "FakeRate", "lnN", SystMap<>::init(1.20));
        
        cb.cp().process({"W"}).channel({"tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                                       "CMS_htt_wjXsec_13TeV", "lnN", SystMap<>::init(1.04));
        
        if(control_region==0){
          cb.cp().process({"W"}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                       "CMS_htt_wjXsec_13TeV", "lnN", SystMap<>::init(1.04));    
        }

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"em","em_2016","em_2017","et","et_2016","et_2017","mt","mt_2016","mt_2017","tt","tt_2016","tt_2017"}).AddSyst(cb,
                                             "CMS_PreFire_13TeV", "shape", SystMap<>::init(1.00));

        // QCD uncerts for em
        // add lnN uncertainties for now        

 
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_em_QCD_BackgroundSubtraction_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_em_QCD_IsoExtrap_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets0_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets0_unc2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets1_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets1_unc2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets2_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1},false).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets2_unc2_13TeV", "shape", SystMap<>::init(1.00));


        //##############################################################################
        //  Fake-Factor Method Uncertainties
        //##############################################################################

//        // FF statistical uncertainties are uncorrelated between all years
//        //
//
//        cb.cp().process({"jetFakes"}).channel({"et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_osss_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_stat_njets0_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_et_qcd_stat_njets0_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_et_qcd_stat_njets1_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_et_qcd_stat_njets1_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_njets2_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_njets2_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_mt_syst_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_mt_syst_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_et_wjets_stat_njets0_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_et_wjets_stat_njets1_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_et_wjets_stat_njets1_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_njets2_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_njets2_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"et_2017","et_2018"}).AddSyst(cb, "ff_et_ttbar_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_ttbar_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_ttbar_stat_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_ttbar_stat_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_sub_syst", "shape", SystMap<>::init(1.00));
//
//
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_osss_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_stat_njets0_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_mt_qcd_stat_njets0_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_mt_qcd_stat_njets1_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_mt_qcd_stat_njets1_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_njets2_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_njets2_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_mt_syst_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_mt_syst_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_mt_wjets_stat_njets0_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_mt_wjets_stat_njets1_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_mt_wjets_stat_njets1_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_njets2_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_njets2_unc2", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_closure_low_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_closure_high_syst", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_stat_unc1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_stat_unc2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_sub_syst", "shape", SystMap<>::init(1.00));
//
//
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_tt_qcd_closure_syst_njets0", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_tt_qcd_closure_syst_njets1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_closure_syst_njets2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2,3,4,5,6},false).AddSyst(cb, "ff_tt_qcd_stat_njets0", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,3,4,5,6},false).AddSyst(cb, "ff_tt_qcd_stat_njets1", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_njets2", "shape", SystMap<>::init(1.00));
//        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_sub_syst", "shape", SystMap<>::init(1.00));
//
//        cb.cp().process({"Wfakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_mc", "lnN", SystMap<>::init(1.3));

        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));

        // only need these uncertainties once using the dR corrected FFs
        cb.cp().process({"jetFakes"}).channel({"tt_2016"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2016", "lnN", SystMap<>::init(1.04));
        cb.cp().process({"jetFakes"}).channel({"tt_2017"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2017", "lnN", SystMap<>::init(1.07));
        cb.cp().process({"jetFakes"}).channel({"tt_2018"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2018", "lnN", SystMap<>::init(1.06));

        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_sub_syst", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1},false).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_dm11", "shape", SystMap<>::init(1.00));

        // put an additional uncertainty on the Wfake comonent which is estimated from MC
        cb.cp().process({"Wfakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_mc", "lnN", SystMap<>::init(1.3));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_sub_syst", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_dm11", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_qcd_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_qcd_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "ff_mt_wjets_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "ff_mt_wjets_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1},false).AddSyst(cb, "ff_mt_wjets_syst_njets2", "shape", SystMap<>::init(1.00));

        // missing template
        // cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_syst", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_sub_syst", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_dm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_dm0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_dm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_dm10", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_dm11", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_qcd_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_qcd_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "ff_et_wjets_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1},false).AddSyst(cb, "ff_et_wjets_syst_njets2", "shape", SystMap<>::init(1.00));

        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_ttbar_syst", "shape", SystMap<>::init(1.00));

        //##############################################################################
        //  DY LO->NLO reweighting, Between no and twice the correction.
        //##############################################################################
        //
        cb.cp().process( {"ZTT","ZJ","ZL","ZLL"}).AddSyst(cb,
                                             "CMS_htt_dyShape_13TeV", "shape", SystMap<>::init(0.1));
        
        
        //##############################################################################
        // Ttbar shape reweighting, Between no and twice the correction
        //##############################################################################
        
        cb.cp().process( {"TTJ","TTT","TT"}).AddSyst(cb,
                                        "CMS_htt_ttbarShape_13TeV", "shape", SystMap<>::init(1.00));
        
        //##############################################################################
        // ZL shape  and electron/muon  to tau fake only in  mt and et channels (updated March 22)
        //##############################################################################

        cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                                         "CMS_ZLShape_et_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                                         "CMS_ZLShape_et_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                         "CMS_ZLShape_mt_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                         "CMS_ZLShape_mt_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
       
        // weighted avarages of recommended tau POG uncertainties provided in bins of eta
       
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).AddSyst(cb, "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.18));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).AddSyst(cb, "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.28));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).AddSyst(cb, "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.23));
        cb.cp().process({"ZL","EWKZ"}).channel({"et_2016"}).AddSyst(cb, "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.24));
        cb.cp().process({"ZL","EWKZ"}).channel({"et_2017"}).AddSyst(cb, "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.21));
        cb.cp().process({"ZL","EWKZ"}).channel({"et_2018"}).AddSyst(cb, "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.15));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2016"}).AddSyst(cb, "CMS_htt_eFakeTau_tt_2016_13TeV", "lnN", SystMap<>::init(1.09));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2017"}).AddSyst(cb, "CMS_htt_eFakeTau_tt_2017_13TeV", "lnN", SystMap<>::init(1.13));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2018"}).AddSyst(cb, "CMS_htt_eFakeTau_tt_2018_13TeV", "lnN", SystMap<>::init(1.11));

 
        //##############################################################################
        // Theoretical Uncertainties on signal (update me)
        //##############################################################################
        // don't use acceptance uncertainties on VBF as there isn't an easy way to get these for the JHU samples (and they are expected to be small for this process)
        // Removed PDF acceptance uncertainties for ggH as these are verysmall compared to PDF uncertainty on XS and scale uncertainty on acceptance/shape
        
        //scale_gg on signal
        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_scale_gg_13TeV", "shape", SystMap<>::init(1.00));
      
        cb.cp().process({"qqH_htt","qqH_htt125"}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_scale_gg_13TeV", "shape", SystMap<>::init(1.00));
 
        //cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).AddSyst(cb,
        //                                     "CMS_FiniteQuarkMass_13TeV", "shape", SystMap<>::init(1.00)); // this uncertainty takes the difference between the finite top-mass dependence and the EFT

        // allow this uncertainty to be unconstrained between the noiminal and +1 sigma (which is difference between top loop and heavy particle loop)
        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).AddSyst(cb,
                                             "CMS_FiniteQuarkMass_13TeV", "shapeU", SystMap<>::init(1.00)); // this uncertainty takes the difference between the finite top-mass dependence and the EFT
        cb.GetParameter("CMS_FiniteQuarkMass_13TeV")->set_range(0, 1);

        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_PS_ISR_ggH_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_PS_FSR_ggH_13TeV", "shape", SystMap<>::init(1.00));

        // shape uncerts for PS for VBF in 2018 where weights exist, lnN for other years (based on 2018 lnN's)
        //cb.cp().process({"qqH_htt125"}).channel({"et_2018","mt_2018","tt_2018","em_2018"}).AddSyst(cb,
        //                                     "CMS_PS_FSR_ggH_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"qqH_htt125"}).channel({"et_2018","mt_2018","tt_2018","em_2018"}).AddSyst(cb,
        //                                     "CMS_PS_ISR_ggH_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.9287)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9938,1.014)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.012,0.9752)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9664,1.017)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.01,0.9745 )); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.012,0.9687)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9797,1.062)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9985,1.006)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.035 )); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.009,1.01  )); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.9857)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.007,0.995 )); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9628,1.023)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9864,1.027)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.9613)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9911,1.096)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.013,0.9898)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.014,0.9532)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9818,1.033)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9975,1.003)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.9959)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.018,0.9712)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.9938)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_FSR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.012,0.9899));
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9988,1.002)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9996,1.001)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9915,1.012)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.9966)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.9986)); 
        cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1,1.002     )); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9971,1.003)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.9979)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1,1.001     )); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9988,1.001)); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998 )); 
        cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.9987)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.001,1     )); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9965,1.005)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9928,1.009)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9842,1.027)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(0.9995,1.001)); 
        cb.cp().channel({"et_2016","et_2017","et_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.9977)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({1}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998 )); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({2}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1,1.001     )); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({3}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.9984)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({4}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.994 )); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({5}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.9961)); 
        cb.cp().channel({"em_2016","em_2017","em_2018"}).bin_id({6}).process({qqH_sig_procs}).AddSyst(cb, "CMS_PS_ISR_VBF_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.9974));


        // extra PS uncerts for VBF to account for the different with the newer samples with Dipole showering
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({1}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.040));
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({2}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.008));
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({3}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.017));
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({4}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.044));
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({5}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.050));
        cb.cp().process({qqH_sig_procs}).channel({"em_2016","em_2017","em_2018"}).bin_id({6}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.022));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({1}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.990));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({2}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.983));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({3}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.998));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({4}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.938));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({5}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.026));
        cb.cp().process({qqH_sig_procs}).channel({"et_2016","et_2017","et_2018"}).bin_id({6}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.021));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.008));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({2}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.014));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({3}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.966));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({4}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.987));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({5}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.047));
        cb.cp().process({qqH_sig_procs}).channel({"mt_2016","mt_2017","mt_2018"}).bin_id({6}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.044));
        cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.127));
        cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({2}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.029));
        cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({3}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(0.984));
        //cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({4}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.000));
        cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({5}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.074));
        cb.cp().process({qqH_sig_procs}).channel({"tt_2016","tt_2017","tt_2018"}).bin_id({6}).AddSyst(cb, "CMS_PS_VBF_withDipoleRecoil_13TeV", "lnN", SystMap<>::init(1.034));

        cb.cp().process(ggH_sig_procs).channel({"et","et_2016","mt","mt_2016","tt","tt_2016","em","em_2016","et_2017","et_2018","mt_2017","mt_2018","tt_2017","tt_2018","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_UE_ggH_13TeV", "shape", SystMap<>::init(1.00));
        
        //    Uncertainty on BR for HTT @ 125 GeV
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_THU", "lnN", SystMap<>::init(1.017));
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_PU_mq", "lnN", SystMap<>::init(1.0099));
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_PU_alphas", "lnN", SystMap<>::init(1.0062));
        
        //    Uncertainty on BR of HWW @ 125 GeV
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_THU", "lnN", SystMap<>::init(1.0099));
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_PU_mq", "lnN", SystMap<>::init(1.0099));
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_PU_alphas", "lnN", SystMap<>::init(1.0066));
        
        
        cb.cp().process(JoinStr({ggH_sig_procs, {"ggH_hww125"}})).AddSyst(cb,"QCDscale_ggH", "lnN", SystMap<>::init(1.039));
        cb.cp().process(JoinStr({qqH_sig_procs, {"qqH_hww125"}})).AddSyst(cb,"QCDscale_qqH", "lnN", SystMap<>::init(1.004));
        cb.cp().process({"wh125_powheg","WH_htt","WH_sm_htt125","WH_ps_htt","WH_mm_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt"}).AddSyst(cb,"QCDscale_qqH", "lnN", SystMap<>::init(1.007));
        cb.cp().process({"wh125_powheg","ZH_htt","ZH_sm_htt125","ZH_ps_htt","ZH_mm_htt","ZH_sm_htt","ZH_ps_htt","ZH_mm_htt"}).AddSyst(cb,"QCDscale_qqH", "lnN", SystMap<>::init(1.038));
        
        cb.cp().process(JoinStr({ggH_sig_procs, {"ggH_hww125"}})).AddSyst(cb,"pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
        cb.cp().process(JoinStr({qqH_sig_procs, {"qqH_hww125"}})).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
        cb.cp().process({"wh125_powheg","WH_htt","WH_sm_htt125","WH_ps_htt","WH_mm_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt"}).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.019));
        cb.cp().process({"zh125_powheg","ZH_htt""ZH_sm_htt125","ZH_ps_htt","ZH_mm_htt","ZH_sm_htt","ZH_ps_htt","ZH_mm_htt"}).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.016));
        
        // jet bin migration uncertainties from: https://arxiv.org/pdf/1610.07922.pdf#subsection.1.4.2.5 (Table 20)
        // For boosted category this is not exclusivly 1 jet events since events with > 1 jets and mjj<300 enter also. So take weighted average of Njets=1 and Njets>=1 uncertainties i.e sigma(boosted) = sigma(njets=1)*(# Njets=1 && boosted)/(# boosted) + sigma(njets>=1)*(#Njets>1 && boosted)/(# boosted)
        // These need to be set properly for MVA approach (placeholders for now)
        
        //cb.cp().AddSyst(cb, "CMS_ggH_mig01", "lnN", SystMap<channel, bin_id, process>::init
        //                ({"em","em_2016","em_2017","em_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
        //                ({"et","et_2016","et_2017","et_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
        //                ({"mt","mt_2016","mt_2017","mt_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
        //                ({"tt","tt_2016","tt_2017","tt_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
        //                
        //                ({"em","em_2016","em_2017","em_2018"},{2},ggH_sig_procs, 1.071)
        //                ({"et","et_2016","et_2017","et_2018"},{2},ggH_sig_procs, 1.071)
        //                ({"mt","mt_2016","mt_2017","mt_2018"},{2},ggH_sig_procs, 1.071)
        //                ({"tt","tt_2016","tt_2017","tt_2018"},{2},ggH_sig_procs, 1.071)
        //                
        //                ({"em","em_2016","em_2017","em_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
        //                ({"et","et_2016","et_2017","et_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
        //                ({"mt","mt_2016","mt_2017","mt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
        //                ({"tt","tt_2016","tt_2017","tt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
        //                );
        //
        //
        //cb.cp().AddSyst(cb, "CMS_ggH_mig12", "lnN", SystMap<channel, bin_id, process>::init 
        //                ({"em","em_2016","em_2017","em_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
        //                ({"et","et_2016","et_2017","et_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
        //                ({"mt","mt_2016","mt_2017","mt_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
        //                ({"tt","tt_2016","tt_2017","tt_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
        //                
        //                ({"em","em_2016","em_2017","em_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
        //                ({"et","et_2016","et_2017","et_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
        //                ({"mt","mt_2016","mt_2017","mt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
        //                ({"tt","tt_2016","tt_2017","tt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
        //                );

        
        //if (ttbar_fit) {
        //    cb.cp().channel({"ttbar_2016","em_2016","et_2016","mt_2016","tt_2016"}).process({"TT","TTT","TTJ"}).AddSyst(cb, "rate_ttbar_2016", "rateParam", SystMap<>::init(1.0));
        //    cb.cp().channel({"ttbar_2017","em_2017","et_2017","mt_2017","tt_2017"}).process({"TT","TTT","TTJ"}).AddSyst(cb, "rate_ttbar_2017", "rateParam", SystMap<>::init(1.0));
        //    cb.cp().channel({"ttbar_2018","em_2018","et_2018","mt_2018","tt_2018"}).process({"TT","TTT","TTJ"}).AddSyst(cb, "rate_ttbar_2018", "rateParam", SystMap<>::init(1.0));
        //    for (auto sys : cb.cp().syst_type({"rateParam"}).syst_name_set()) {
        //        if(sys.find("rate_ttbar") != std::string::npos) cb.GetParameter(sys)->set_range(0.8, 1.2); 
        //    }
        //}
        
    }
}
