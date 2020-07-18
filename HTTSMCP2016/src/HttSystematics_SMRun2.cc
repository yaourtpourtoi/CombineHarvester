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
        
        
        std::vector<std::string> sig_procs = {"ggH_htt","qqH_htt","WH_htt","ZH_htt","ggHsm_htt", "ggHps_htt", "ggHmm_htt","qqHsm_htt", "qqHps_htt", "qqHmm_htt","qqH_htt125","qqHsm_htt125", "qqHps_htt", "qqHmm_htt","WH_htt125","ZH_htt125","WHsm_htt125","ZHsm_htt125", "WHps_htt","ZHps_htt","WHmm_htt","ZHmm_htt","WHsm_htt","ZHsm_htt","WHps_htt","ZHps_htt","WHmm_htt","ZHmm_htt", "ggHsm_jhu_htt","ggHps_jhu_htt","ggHmm_jhu_htt","ggH_ph_htt"};
        std::vector<std::string> ggH_sig_procs = {"ggH_htt","ggHsm_htt", "ggHps_htt", "ggHmm_htt","ggHsm_jhu_htt","ggHps_jhu_htt","ggHmm_jhu_htt","ggH_ph_htt"};
        std::vector<std::string> qqH_sig_procs = {"qqH_htt","qqHsm_htt", "qqHps_htt", "qqHmm_htt", "qqH_htt125","qqHsm_htt125", "qqHps_htt", "qqHmm_htt"};
        
        // N.B. when adding this list of backgrounds to a nuisance, only
        // the backgrounds that are included in the background process
        // defined in MorphingSM2016.cpp are included in the actual DCs
        // This is a list of all MC based backgrounds
        // QCD is explicitly excluded
        std::vector<std::string> all_mc_bkgs = {
            "ZL","ZLL","ZJ","ZTT","TTJ","TTT","TT",
            "W","VV","VVT","VVJ",
            "ggH_hww125","qqH_hww125","EWKZ"};
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

        // uncorrelated parts 2.2%, 2.0%, 1.5%
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_2016_13TeV", "lnN", SystMap<>::init(1.022));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_2017_13TeV", "lnN", SystMap<>::init(1.020));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_2018_13TeV", "lnN", SystMap<>::init(1.015));

        // correlated parts, 1.2%, 1.1%, 2.0%
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","et","et_2016","mt","mt_2016","em","em_2016","ttbar","ttbar_2016"}).AddSyst(cb,
                                            "lumi_13TeV", "lnN", SystMap<>::init(1.012));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2017","mt_2017","et_2017","em_2017","ttbar_2017"}).AddSyst(cb,
                                            "lumi_13TeV", "lnN", SystMap<>::init(1.011));
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt_2018","mt_2018","et_2018","em_2018","ttbar_2018"}).AddSyst(cb,
                                            "lumi_13TeV", "lnN", SystMap<>::init(1.020));

        //##############################################################################
        //  trigger   
        //##############################################################################
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_mt_13TeV", "lnN", SystMap<>::init(1.02));
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_et_13TeV", "lnN", SystMap<>::init(1.02)); 
        
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs,embed})).channel({"em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                             "CMS_eff_trigger_em_13TeV", "lnN", SystMap<>::init(1.02)); 


        // use lN uncertainty for tt trigger until the differential ones are fixed

//        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
//                                             "CMS_eff_trigger_tt_13TeV", "lnN", SystMap<>::init(1.1));

        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_trg_DM11_13TeV", "shape", SystMap<>::init(1.00));


        // additional uncertainties due to tau SF on cross-triggers - decorrelated for embedded samples in Morphing code
        //  Add back when this is working properly (i.e when you have fixed issue with nan values in weights)
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM0_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM1_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_mt_DM11_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM0_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM1_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_DM11_13TeV", "shape", SystMap<>::init(1.00));        

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


//CMS_eff_m, CMS_eff_e, CMS_eff_t_mt_13TeV, CMS_eff_t_et_13TeV, CMS_eff_t_tt_13TeV

        // due to different treatments of embedding and MC uncertainties for the tau ID they are included seperatly for now

        cb.cp().process(embed).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_13TeV", "lnN", SystMap<>::init(1.05));

        cb.cp().process(embed).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_13TeV", "lnN", SystMap<>::init(1.1));

        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM10_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_DM11_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin3_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin4_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,real_tau_mc_bkgs})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,"CMS_eff_t_bin5_13TeV", "shape", SystMap<>::init(1.00));

        // the uncorrelated part by channel is due to anti-electron and anti-muon disctiminant, it is the same for MC and embedding so decorrelated in the Morphing       
 
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_mt_13TeV", "lnN", SystMap<>::init(1.01));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                             "CMS_eff_t_et_13TeV", "lnN", SystMap<>::init(1.01));

        // 2 real taus
        cb.cp().process(JoinStr({sig_procs, {"ZTT","VVT","TTT","EWKZ"}, embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.02));

        // 1 real tau
        cb.cp().process(JoinStr({{"TTJ","ZJ","VVJ","W","Wfakes","TTfakes"}})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.01));


//CMS_eff_t_DM_13TeV
        
        //##############################################################################
        //  b tag and mistag rate  efficiencies 
        //##############################################################################
 
        // update number for 2017 
        // real uncerts for TT and VV only (others are small) 
        cb.cp().AddSyst(cb,
          "CMS_eff_b_13TeV", "lnN", SystMapAsymm<channel,bin_id,process>::init
          ({"et","et_2016"}, {1}, {"TTT"}, 1.020,0.969)
          ({"et","et_2016"}, {2}, {"TTT"}, 1.034,0.968)
          ({"et","et_2016"}, {3}, {"TTT"}, 1.025,0.932)
          ({"et","et_2016"}, {4}, {"TTT"}, 1.008,0.969)
          ({"et","et_2016"}, {5}, {"TTT"}, 1.032,0.963)
          ({"et","et_2016"}, {6}, {"TTT"}, 1.015,0.942)
          ({"et","et_2016"}, {1}, {"VVT"}, 1.001,1.0)
          ({"et","et_2016"}, {2}, {"VVT"}, 1.010,0.989)
          ({"et","et_2016"}, {3}, {"VVT"}, 1.024,0.987)
          ({"et","et_2016"}, {4}, {"VVT"}, 1.000,0.983)
          //({"et","et_2016"}, {5}, {"VVT"}, 1.000,1.000)
          //({"et","et_2016"}, {6}, {"VVT"}, 1.000,1.000)
          ({"mt","mt_2016"}, {1}, {"TTT"}, 1.024,0.975)
          ({"mt","mt_2016"}, {2}, {"TTT"}, 1.030,0.967)
          ({"mt","mt_2016"}, {3}, {"TTT"}, 1.069,0.967)
          ({"mt","mt_2016"}, {4}, {"TTT"}, 1.040,0.964)
          ({"mt","mt_2016"}, {5}, {"TTT"}, 1.028,0.968)
          ({"mt","mt_2016"}, {6}, {"TTT"}, 1.060,0.952)
          ({"mt","mt_2016"}, {1}, {"VVT"}, 1.0001,0.999) 
          ({"mt","mt_2016"}, {2}, {"VVT"}, 1.010,0.994)
          ({"mt","mt_2016"}, {3}, {"VVT"}, 1.015,1.000)
          ({"mt","mt_2016"}, {4}, {"VVT"}, 1.014,0.987)
          ({"mt","mt_2016"}, {5}, {"VVT"}, 1.005,1.000)
          ({"mt","mt_2016"}, {6}, {"VVT"}, 1.024,0.986)
          ({"em","em_2016"}, {1}, {"TT"}, 1.016,0.982)
          ({"em","em_2016"}, {2}, {"TT"}, 1.030,0.968)
          ({"em","em_2016"}, {3}, {"TT"}, 1.033,0.968)
          ({"em","em_2016"}, {4}, {"TT"}, 1.032,0.957)
          ({"em","em_2016"}, {5}, {"TT"}, 1.032,0.975)
          ({"em","em_2016"}, {6}, {"TT"}, 1.032,0.969)
          ({"em","em_2016"}, {1}, {"VV"}, 1.001, 0.999)
          ({"em","em_2016"}, {2}, {"VV"}, 1.008,0.992)
          ({"em","em_2016"}, {3}, {"VV"}, 1.013,0.986)
          ({"em","em_2016"}, {4}, {"VV"}, 1.008,0.988)
          ({"em","em_2016"}, {5}, {"VV"}, 1.010,0.992)
          ({"em","em_2016"}, {6}, {"VV"}, 1.011,0.988)
          ({"em_2017"}, {1}, {"TT"}, 1.032, 0.973)
          ({"em_2017"}, {1}, {"VV"}, 1.002, 0.998)
          ({"em_2017"}, {4}, {"TT"}, 1.042, 0.963)
          ({"em_2017"}, {4}, {"VV"}, 1.018, 0.980)
          ({"em_2017"}, {3}, {"TT"}, 1.036, 0.962)
          ({"em_2017"}, {3}, {"VV"}, 1.014, 0.985)
          ({"em_2017"}, {6}, {"TT"}, 1.041, 0.959)
          ({"em_2017"}, {6}, {"VV"}, 1.017, 0.979)
          ({"em_2017"}, {5}, {"TT"}, 1.050, 0.965)
          ({"em_2017"}, {5}, {"VV"}, 1.007, 0.980)
          ({"em_2017"}, {2}, {"TT"}, 1.036, 0.963)
          ({"em_2017"}, {2}, {"VV"}, 1.011, 0.991)
          ({"et_2017"}, {1}, {"TTT"}, 1.049, 0.981)
          ({"et_2017"}, {1}, {"VVT"}, 1.001, 0.998)
          ({"et_2017"}, {4}, {"TTT"}, 1.032, 0.975)
          ({"et_2017"}, {4}, {"VVT"}, 1.000, 0.974)
          ({"et_2017"}, {3}, {"TTT"}, 1.052, 0.954)
          ({"et_2017"}, {3}, {"VVT"}, 1.037, 0.964)
          ({"et_2017"}, {6}, {"TTT"}, 1.015, 0.925)
          ({"et_2017"}, {6}, {"VVT"}, 1.026, 1.002)
          ({"et_2017"}, {5}, {"TTT"}, 1.020, 0.977)
          //({"et_2017"}, {5}, {"VVT"}, 1.000, 1.000)
          ({"et_2017"}, {2}, {"TTT"}, 1.047, 0.957)
          ({"et_2017"}, {2}, {"VVT"}, 1.010, 0.988)
          ({"mt_2017"}, {1}, {"TTT"}, 1.005, 0.972)
          ({"mt_2017"}, {1}, {"VVT"}, 1.002, 0.998)
          ({"mt_2017"}, {4}, {"TTT"}, 1.047, 0.968)
          ({"mt_2017"}, {4}, {"VVT"}, 1.040, 0.979)
          ({"mt_2017"}, {3}, {"TTT"}, 1.042, 0.944)
          ({"mt_2017"}, {3}, {"VVT"}, 1.008, 0.974)
          ({"mt_2017"}, {6}, {"TTT"}, 1.074, 0.964)
          ({"mt_2017"}, {6}, {"VVT"}, 1.058, 0.973)
          ({"mt_2017"}, {5}, {"TTT"}, 1.057, 0.981)
          ({"mt_2017"}, {5}, {"VVT"}, 1.003, 0.961)
          ({"mt_2017"}, {2}, {"TTT"}, 1.038, 0.956)
          ({"mt_2017"}, {2}, {"VVT"}, 1.009, 0.988)

          ({"em_2018"}, {1}, {"TT"}, 1.032, 0.973)
          ({"em_2018"}, {1}, {"VV"}, 1.002, 0.998)
          ({"em_2018"}, {4}, {"TT"}, 1.042, 0.963)
          ({"em_2018"}, {4}, {"VV"}, 1.018, 0.980)
          ({"em_2018"}, {3}, {"TT"}, 1.036, 0.962)
          ({"em_2018"}, {3}, {"VV"}, 1.014, 0.985)
          ({"em_2018"}, {6}, {"TT"}, 1.041, 0.959)
          ({"em_2018"}, {6}, {"VV"}, 1.017, 0.979)
          ({"em_2018"}, {5}, {"TT"}, 1.050, 0.965)
          ({"em_2018"}, {5}, {"VV"}, 1.007, 0.980)
          ({"em_2018"}, {2}, {"TT"}, 1.036, 0.963)
          ({"em_2018"}, {2}, {"VV"}, 1.011, 0.991)
          ({"et_2018"}, {1}, {"TTT"}, 1.049, 0.981)
          ({"et_2018"}, {1}, {"VVT"}, 1.001, 0.998)
          ({"et_2018"}, {4}, {"TTT"}, 1.032, 0.975)
          ({"et_2018"}, {4}, {"VVT"}, 1.000, 0.974)
          ({"et_2018"}, {3}, {"TTT"}, 1.052, 0.954)
          ({"et_2018"}, {3}, {"VVT"}, 1.037, 0.964)
          ({"et_2018"}, {6}, {"TTT"}, 1.015, 0.925)
          ({"et_2018"}, {6}, {"VVT"}, 1.026, 1.002)
          ({"et_2018"}, {5}, {"TTT"}, 1.020, 0.977)
          //({"et_2018"}, {5}, {"VVT"}, 1.000, 1.000)
          ({"et_2018"}, {2}, {"TTT"}, 1.047, 0.957)
          ({"et_2018"}, {2}, {"VVT"}, 1.010, 0.988)
          ({"mt_2018"}, {1}, {"TTT"}, 1.005, 0.972)
          ({"mt_2018"}, {1}, {"VVT"}, 1.002, 0.998)
          ({"mt_2018"}, {4}, {"TTT"}, 1.047, 0.968)
          ({"mt_2018"}, {4}, {"VVT"}, 1.040, 0.979)
          ({"mt_2018"}, {3}, {"TTT"}, 1.042, 0.944)
          ({"mt_2018"}, {3}, {"VVT"}, 1.008, 0.974)
          ({"mt_2018"}, {6}, {"TTT"}, 1.074, 0.964)
          ({"mt_2018"}, {6}, {"VVT"}, 1.058, 0.973)
          ({"mt_2018"}, {5}, {"TTT"}, 1.057, 0.981)
          ({"mt_2018"}, {5}, {"VVT"}, 1.003, 0.961)
          ({"mt_2018"}, {2}, {"TTT"}, 1.038, 0.956)
          ({"mt_2018"}, {2}, {"VVT"}, 1.009, 0.988)

          // the commented numbers are for deepCSV (above numbers for CSVv2)
          //({"mt_2017"}, {1}, {"TTT"}, 1.026, 0.982)
          //({"mt_2017"}, {1}, {"VVT"}, 1.002, 0.999)
          //({"mt_2017"}, {4}, {"TTT"}, 1.021, 0.965)
          //({"mt_2017"}, {4}, {"VVT"}, 1.039, 0.986)
          //({"mt_2017"}, {3}, {"TTT"}, 1.072, 0.971)
          //({"mt_2017"}, {3}, {"VVT"}, 1.022, 0.982)
          //({"mt_2017"}, {6}, {"TTT"}, 1.072, 0.972)
          //({"mt_2017"}, {6}, {"VVT"}, 1.028, 0.963)
          //({"mt_2017"}, {5}, {"TTT"}, 1.047, 0.989)
          //({"mt_2017"}, {5}, {"VVT"}, 1.023, 0.955)
          //({"mt_2017"}, {2}, {"TTT"}, 1.052, 0.950)
          //({"mt_2017"}, {2}, {"VVT"}, 1.012, 0.986)
          //({"et_2017"}, {1}, {"TTT"}, 1.009, 0.976)
          //({"et_2017"}, {1}, {"VVT"}, 1.001, 0.999)
          //({"et_2017"}, {4}, {"TTT"}, 1.065, 0.968)
          //({"et_2017"}, {4}, {"VVT"}, 1.019, 0.975)
          //({"et_2017"}, {3}, {"TTT"}, 1.062, 0.906)
          //({"et_2017"}, {3}, {"VVT"}, 1.043, 1.000)
          //({"et_2017"}, {6}, {"TTT"}, 1.068, 0.947)
          //({"et_2017"}, {6}, {"VVT"}, 1.024, 1.001)
          //({"et_2017"}, {5}, {"TTT"}, 1.026, 0.988)
          //({"et_2017"}, {5}, {"VVT"}, 1.000, 1.000)
          //({"et_2017"}, {2}, {"TTT"}, 1.046, 0.950)
          //({"et_2017"}, {2}, {"VVT"}, 1.016, 0.987)
          //({"em_2017"}, {1}, {"TT"}, 1.029, 0.980)
          //({"em_2017"}, {1}, {"VV"}, 1.002, 0.999)
          //({"em_2017"}, {4}, {"TT"}, 1.055, 0.948)
          //({"em_2017"}, {4}, {"VV"}, 1.021, 0.976)
          //({"em_2017"}, {3}, {"TT"}, 1.047, 0.948)
          //({"em_2017"}, {3}, {"VV"}, 1.020, 0.987)
          //({"em_2017"}, {6}, {"TT"}, 1.050, 0.946)
          //({"em_2017"}, {6}, {"VV"}, 1.022, 0.982)
          //({"em_2017"}, {5}, {"TT"}, 1.050, 0.951)
          //({"em_2017"}, {5}, {"VV"}, 1.017, 0.979)
          //({"em_2017"}, {2}, {"TT"}, 1.044, 0.957)
          //({"em_2017"}, {2}, {"VV"}, 1.010, 0.989)
          //
          //({"mt_2018"}, {1}, {"TTT"}, 1.026, 0.982)
          //({"mt_2018"}, {1}, {"VVT"}, 1.002, 0.999)
          //({"mt_2018"}, {4}, {"TTT"}, 1.021, 0.965)
          //({"mt_2018"}, {4}, {"VVT"}, 1.039, 0.986)
          //({"mt_2018"}, {3}, {"TTT"}, 1.072, 0.971)
          //({"mt_2018"}, {3}, {"VVT"}, 1.022, 0.982)
          //({"mt_2018"}, {6}, {"TTT"}, 1.072, 0.972)
          //({"mt_2018"}, {6}, {"VVT"}, 1.028, 0.963)
          //({"mt_2018"}, {5}, {"TTT"}, 1.047, 0.989)
          //({"mt_2018"}, {5}, {"VVT"}, 1.023, 0.955)
          //({"mt_2018"}, {2}, {"TTT"}, 1.052, 0.950)
          //({"mt_2018"}, {2}, {"VVT"}, 1.012, 0.986)
          //({"et_2018"}, {1}, {"TTT"}, 1.009, 0.976)
          //({"et_2018"}, {1}, {"VVT"}, 1.001, 0.999)
          //({"et_2018"}, {4}, {"TTT"}, 1.065, 0.968)
          //({"et_2018"}, {4}, {"VVT"}, 1.019, 0.975)
          //({"et_2018"}, {3}, {"TTT"}, 1.062, 0.906)
          //({"et_2018"}, {3}, {"VVT"}, 1.043, 1.000)
          //({"et_2018"}, {6}, {"TTT"}, 1.068, 0.947)
          //({"et_2018"}, {6}, {"VVT"}, 1.024, 1.001)
          //({"et_2018"}, {5}, {"TTT"}, 1.026, 0.988)
          //({"et_2018"}, {5}, {"VVT"}, 1.000, 1.000)
          //({"et_2018"}, {2}, {"TTT"}, 1.046, 0.950)
          //({"et_2018"}, {2}, {"VVT"}, 1.016, 0.987)
          //({"em_2018"}, {1}, {"TT"}, 1.029, 0.980)
          //({"em_2018"}, {1}, {"VV"}, 1.002, 0.999)
          //({"em_2018"}, {4}, {"TT"}, 1.055, 0.948)
          //({"em_2018"}, {4}, {"VV"}, 1.021, 0.976)
          //({"em_2018"}, {3}, {"TT"}, 1.047, 0.948)
          //({"em_2018"}, {3}, {"VV"}, 1.020, 0.987)
          //({"em_2018"}, {6}, {"TT"}, 1.050, 0.946)
          //({"em_2018"}, {6}, {"VV"}, 1.022, 0.982)
          //({"em_2018"}, {5}, {"TT"}, 1.050, 0.951)
          //({"em_2018"}, {5}, {"VV"}, 1.017, 0.979)
          //({"em_2018"}, {2}, {"TT"}, 1.044, 0.957)
          //({"em_2018"}, {2}, {"VV"}, 1.010, 0.989)
          
          ({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}, {100,101,102,103,104,107,108}, {"TTT"}, 1.05, 0.95)
          ({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}, {105,106}, {"TTT"}, 1.06, 0.91)
          ({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}, {100,101,102,103,104,107,108}, {"VVT"}, 1.02, 0.98)
          ({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}, {105,106}, {"VVT"}, 1.04, 0.96)
        );

        cb.cp().process({"TTT","TT","VVT","VV"}).bin_id({31,32,33,34,35,36,37,41,42,43,44,45,46,47}).AddSyst(cb,
                                             "CMS_eff_b_13TeV", "shape", SystMap<>::init(1.00));
        
        //##############################################################################
        //  Electron, muon and tau energy Scale
        //##############################################################################
        
        // Add back later!
        //cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed/*, {"jetFakes", "QCD"}*/})).channel({"et","et_2016","et_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
        //                                     "CMS_scale_e_13TeV", "shape", SystMap<>::init(1.00));
       

 
        // Decay Mode based TES Settings
        //cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
        //                                        "CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
        //                                        "CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed, /*{"jetFakes"}*/})).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
        //                                        "CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));

        // Muon 
        // ES Add back after!
        //cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed, /*{"jetFakes", "QCD"}*/})).channel({"mt","mt_2016","mt_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
        //                                     "CMS_scale_mu_13TeV", "shape", SystMap<>::init(1.00));

        //##############################################################################
        //  Embedded uncertainty on ttbar contamination (and VV contamination)
        //##############################################################################        
        cb.cp().process(embed).AddSyst(cb,"CMS_ttbar_embeded_13TeV", "shape", SystMap<>::init(1.00));
 
        //##############################################################################
        //  jet and met energy Scale
        //##############################################################################


        // MET Systematic shapes - recoil uncertainties for recoil corrected met, unclustered energy uncertainty for samples with no recoil correction, jes uncertainties propogated to met for samples with no recoil correction
        /*cb.cp().process({"TT","TTJ","TTT","VV","VVJ","VVT"}).channel({"et","et_2016","mt","mt_2016","tt","tt_2016"},false).AddSyst(cb,
                                                  "CMS_scale_met_unclustered_13TeV", "shape", SystMap<>::init(1.00));


        cb.cp().process(JoinStr({sig_procs, {"ZTT","ZLL","ZL","ZJ","EWKZ","W"}})).AddSyst(cb,
                                                  "CMS_htt_boson_reso_met_13TeV", "shape", SystMap<>::init(1.00)); 
        cb.cp().process(JoinStr({sig_procs, {"ZTT","ZLL","ZL","ZJ","EWKZ","W"}})).AddSyst(cb,
                                                  "CMS_htt_boson_scale_met_13TeV", "shape", SystMap<>::init(1.00));      */
 

      

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).AddSyst(cb,"CMS_scale_j_13TeV", "shape", SystMap<>::init(1.00));
 
        
        //##############################################################################
        //  Background normalization uncertainties
        //##############################################################################
        
        //   Diboson  Normalisation - fully correlated
        cb.cp().process({"VV","VVT","VVJ"}).AddSyst(cb,
                                        "CMS_htt_vvXsec_13TeV", "lnN", SystMap<>::init(1.05));

        cb.cp().process({"ZTT","ZJ","ZL","ZLL"}).AddSyst(cb,
                                        "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.04));        
 
        cb.cp().process({"EWKZ"}).AddSyst(cb,
                                        "CMS_htt_ewkzXsec_13TeV", "lnN", SystMap<>::init(1.04));

        //if (! ttbar_fit){
        //   ttbar Normalisation - fully correlated
	    cb.cp().process({"TT","TTT","TTJ"}).AddSyst(cb,
					  "CMS_htt_tjXsec_13TeV", "lnN", SystMap<>::init(1.06));
        //}

        // W norm, just for em, tt and the mm region where MC norm is from MC
        
        cb.cp().process({"W"}).channel({"em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                                       "CMS_htt_jetFakeLep_13TeV", "lnN", SystMap<>::init(1.20));
        
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
                                                       "CMS_$CHANNEL_QCD_$BIN_13TeV", "lnN", SystMap<>::init(1.15));

 
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_em_QCD_BackgroundSubtraction_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_em_QCD_IsoExtrap_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets0_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({1}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets0_unc2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets1_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets1_unc2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2,3,4,5,6}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets2_unc1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process({"QCD"}).channel({"em","em_2016","em_2017","em_2018"}).bin_id({2,3,4,5,6}).AddSyst(cb,
                                             "CMS_em_QCD_stat_njets2_unc2_13TeV", "shape", SystMap<>::init(1.00));


        //##############################################################################
        //  Fake-Factor Method Uncertainties
        //##############################################################################

        // FF statistical uncertainties are uncorrelated between all years
        //
        // add lnN uncertainties fow now!
        //
        cb.cp().process({"jetFakes"}).AddSyst(cb,
                                                       "CMS_$CHANNEL_FF_$BIN_13TeV", "lnN", SystMap<>::init(1.15));

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

       // cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
       //                                                  "CMS_ZLShape_et_1prong_13TeV", "shape", SystMap<>::init(1.00));
       // cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
       //                                                  "CMS_ZLShape_et_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));

       // cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
       //                                                  "CMS_ZLShape_mt_1prong_13TeV", "shape", SystMap<>::init(1.00));
       // cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
       //                                                  "CMS_ZLShape_mt_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
       
        // weighted avarages of recommended tau POG uncertainties provided in bins of eta
        cb.cp().process({"ZL","EWKZ"}).channel({"mt","mt_2016"}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.27));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.27));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.27));
        cb.cp().process({"ZL","EWKZ"}).channel({"et","et_2016"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt","tt_2016"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"et_2017"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"et_2018"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2017"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.23));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2018"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.23));
        
        //##############################################################################
        // Theoretical Uncertainties on signal (update me)
        //##############################################################################
        // don't use acceptance uncertainties on VBF as there isn't an easy way to get these for the JHU samples (and they are expected to be small for this process)
        // Removed PDF acceptance uncertainties for ggH as these are verysmall compared to PDF uncertainty on XS and scale uncertainty on acceptance/shape
        
        //scale_gg on signal
        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).bin_id({2,3,4,5,6,31,32,33,34,35,36,37,41,42,43,44,45,46,47,48,49}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_scale_gg_13TeV", "shape", SystMap<>::init(1.00));
       
        cb.cp().AddSyst(cb,
          "CMS_scale_gg_13TeV", "lnN", SystMapAsymm<channel,bin_id,process>::init
          ({"tt","tt_2016"}, {1}, {"ggHps_htt"}, 0.937, 1.063)
          ({"tt","tt_2016"}, {1}, {"ggHsm_htt"}, 0.928, 1.069)
          ({"tt","tt_2016"}, {1}, {"ggHmm_htt"}, 0.942, 1.059)
          ({"mt","mt_2016"}, {1}, {"ggHps_htt"}, 0.948, 1.054)
          ({"mt","mt_2016"}, {1}, {"ggHsm_htt"}, 0.939, 1.061)
          ({"mt","mt_2016"}, {1}, {"ggHmm_htt"}, 0.946, 1.056)
          ({"et","et_2016"}, {1}, {"ggHps_htt"}, 0.946, 1.056)
          ({"et","et_2016"}, {1}, {"ggHsm_htt"}, 0.966, 1.042)
          ({"et","et_2016"}, {1}, {"ggHmm_htt"}, 0.929, 1.069)
          ({"em","em_2016"}, {1}, {"ggHps_htt"}, 0.947, 1.056)
          ({"em","em_2016"}, {1}, {"ggHsm_htt"}, 0.942, 1.058)
          ({"em","em_2016"}, {1}, {"ggHmm_htt"}, 0.950, 1.053)
          ({"em_2017"}, {1}, {"ggHps_htt"}, 0.975, 1.027)
          ({"em_2017"}, {1}, {"ggHmm_htt"}, 0.963, 1.037)
          ({"em_2017"}, {1}, {"ggHsm_htt"}, 0.954, 1.044)
          ({"et_2017"}, {1}, {"ggHps_htt"}, 0.963, 1.036)
          ({"et_2017"}, {1}, {"ggHmm_htt"}, 0.973, 1.031)
          ({"et_2017"}, {1}, {"ggHsm_htt"}, 0.946, 1.050)
          ({"mt_2017"}, {1}, {"ggHps_htt"}, 0.987, 1.019)
          ({"mt_2017"}, {1}, {"ggHmm_htt"}, 0.966, 1.035)
          ({"mt_2017"}, {1}, {"ggHsm_htt"}, 0.952, 1.047)
          ({"tt_2017"}, {1}, {"ggHps_htt"}, 0.975, 1.029)
          ({"tt_2017"}, {1}, {"ggHmm_htt"}, 0.966, 1.035)
          ({"tt_2017"}, {1}, {"ggHsm_htt"}, 0.956, 1.043)

          ({"em_2018"}, {1}, {"ggHps_htt"}, 0.975, 1.027)
          ({"em_2018"}, {1}, {"ggHmm_htt"}, 0.963, 1.037)
          ({"em_2018"}, {1}, {"ggHsm_htt"}, 0.954, 1.044)
          ({"et_2018"}, {1}, {"ggHps_htt"}, 0.963, 1.036)
          ({"et_2018"}, {1}, {"ggHmm_htt"}, 0.973, 1.031)
          ({"et_2018"}, {1}, {"ggHsm_htt"}, 0.946, 1.050)
          ({"mt_2018"}, {1}, {"ggHps_htt"}, 0.987, 1.019)
          ({"mt_2018"}, {1}, {"ggHmm_htt"}, 0.966, 1.035)
          ({"mt_2018"}, {1}, {"ggHsm_htt"}, 0.952, 1.047)
          ({"tt_2018"}, {1}, {"ggHps_htt"}, 0.975, 1.029)
          ({"tt_2018"}, {1}, {"ggHmm_htt"}, 0.966, 1.035)
          ({"tt_2018"}, {1}, {"ggHsm_htt"}, 0.956, 1.043)
 
        ); 
       
        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).AddSyst(cb,
                                             "CMS_FiniteQuarkMass_13TeV", "shape", SystMap<>::init(1.00)); // this uncertainty takes the difference between the finite top-mass dependence and the EFT


        // PS uncertainty affects njets and pT distribution so is included as a shape uncertainty for the boosted category where pT is fitted
        cb.cp().process(ggH_sig_procs).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_PS_ggH_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(ggH_sig_procs).AddSyst(cb,
                                             "CMS_UE_ggH_13TeV", "shape", SystMap<>::init(1.00));
        
        //    Uncertainty on BR for HTT @ 125 GeV
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_THU", "lnN", SystMap<>::init(1.017));
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_PU_mq", "lnN", SystMap<>::init(1.0099));
        cb.cp().process(sig_procs).AddSyst(cb,"BR_htt_PU_alphas", "lnN", SystMap<>::init(1.0062));
        
        //    Uncertainty on BR of HWW @ 125 GeV
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_THU", "lnN", SystMap<>::init(1.0099));
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_PU_mq", "lnN", SystMap<>::init(1.0099));
        cb.cp().process({"ggH_hww125","qqH_hww125"}).AddSyst(cb,"BR_hww_PU_alphas", "lnN", SystMap<>::init(1.0066));
        
        
        cb.cp().process(JoinStr({ggH_sig_procs, {"ggH_hww125"}})).AddSyst(cb,"QCDScale_ggH", "lnN", SystMap<>::init(1.039));
        cb.cp().process(JoinStr({qqH_sig_procs, {"qqH_hww125"}})).AddSyst(cb,"QCDScale_qqH", "lnN", SystMap<>::init(1.004));
        cb.cp().process({"WH_htt125","WH_htt","WHsm_htt125","WHps_htt","WHmm_htt","WHsm_htt","WHps_htt","WHmm_htt"}).AddSyst(cb,"QCDScale_WH", "lnN", SystMap<>::init(1.007));
        cb.cp().process({"WH_htt125","ZH_htt","ZHsm_htt125","ZHps_htt","ZHmm_htt","ZHsm_htt","ZHps_htt","ZHmm_htt"}).AddSyst(cb,"QCDScale_ZH", "lnN", SystMap<>::init(1.038));
        
        cb.cp().process(JoinStr({ggH_sig_procs, {"ggH_hww125"}})).AddSyst(cb,"pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
        cb.cp().process(JoinStr({qqH_sig_procs, {"qqH_hww125"}})).AddSyst(cb,"pdf_Higgs_qq", "lnN", SystMap<>::init(1.021));
        cb.cp().process({"WH_htt125","WH_htt","WHsm_htt125","WHps_htt","WHmm_htt","WHsm_htt","WHps_htt","WHmm_htt"}).AddSyst(cb,"pdf_Higgs_WH", "lnN", SystMap<>::init(1.019));
        cb.cp().process({"ZH_htt125","ZH_htt""ZHsm_htt125","ZHps_htt","ZHmm_htt","ZHsm_htt","ZHps_htt","ZHmm_htt"}).AddSyst(cb,"pdf_Higgs_ZH", "lnN", SystMap<>::init(1.016));
        
        // jet bin migration uncertainties from: https://arxiv.org/pdf/1610.07922.pdf#subsection.1.4.2.5 (Table 20)
        // For boosted category this is not exclusivly 1 jet events since events with > 1 jets and mjj<300 enter also. So take weighted average of Njets=1 and Njets>=1 uncertainties i.e sigma(boosted) = sigma(njets=1)*(# Njets=1 && boosted)/(# boosted) + sigma(njets>=1)*(#Njets>1 && boosted)/(# boosted)
        // These need to be set properly for MVA approach (placeholders for now)
        
        cb.cp().AddSyst(cb, "CMS_ggH_mig01", "lnN", SystMap<channel, bin_id, process>::init
                        ({"em","em_2016","em_2017","em_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
                        ({"et","et_2016","et_2017","et_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
                        ({"mt","mt_2016","mt_2017","mt_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
                        ({"tt","tt_2016","tt_2017","tt_2018"},{1,31,32,33,34,35,36,37},ggH_sig_procs, 0.959)
                        
                        ({"em","em_2016","em_2017","em_2018"},{2},ggH_sig_procs, 1.071)
                        ({"et","et_2016","et_2017","et_2018"},{2},ggH_sig_procs, 1.071)
                        ({"mt","mt_2016","mt_2017","mt_2018"},{2},ggH_sig_procs, 1.071)
                        ({"tt","tt_2016","tt_2017","tt_2018"},{2},ggH_sig_procs, 1.071)
                        
                        ({"em","em_2016","em_2017","em_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
                        ({"et","et_2016","et_2017","et_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
                        ({"mt","mt_2016","mt_2017","mt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
                        ({"tt","tt_2016","tt_2017","tt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.036)
                        );
        
        
        cb.cp().AddSyst(cb, "CMS_ggH_mig12", "lnN", SystMap<channel, bin_id, process>::init 
                        ({"em","em_2016","em_2017","em_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
                        ({"et","et_2016","et_2017","et_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
                        ({"mt","mt_2016","mt_2017","mt_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
                        ({"tt","tt_2016","tt_2017","tt_2018"},{2,31,32,33,34,35,36,37},ggH_sig_procs, 0.986)
                        
                        ({"em","em_2016","em_2017","em_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
                        ({"et","et_2016","et_2017","et_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
                        ({"mt","mt_2016","mt_2017","mt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
                        ({"tt","tt_2016","tt_2017","tt_2018"},{3,4,5,6,41,42,43,44,45,46,47,48,49},ggH_sig_procs, 1.145)
                        );

        
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
