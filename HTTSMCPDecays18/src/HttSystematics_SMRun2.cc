#include "CombineHarvester/HTTSMCPDecays18/interface/HttSystematics_SMRun2.h"
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
        
        
        std::vector<std::string> sig_procs = {
            "ggH_htt","qqH_htt","WH_htt","ZH_htt","ggH_sm_htt",
            "ggH_ps_htt", "ggH_mm_htt","qqH_sm_htt", "qqH_ps_htt",
            "qqH_mm_htt","qqH_htt125","qqH_sm_htt", "qqH_ps_htt",
            "qqH_mm_htt","WH_htt125","ZH_htt125","WH_sm_htt","ZH_sm_htt",
            "WH_ps_htt","ZH_ps_htt","WH_mm_htt","ZH_mm_htt","WH_sm_htt",
            "ZH_sm_htt","WH_ps_htt","ZH_ps_htt","WH_mm_htt","ZH_mm_htt"};
        std::vector<std::string> ggH_sig_procs = {"ggH_htt","ggH_sm_htt", "ggH_ps_htt",
            "ggH_mm_htt"};
        std::vector<std::string> qqH_sig_procs = {"qqH_htt","qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt",
            "qqH_htt125","qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt"};
        std::vector<std::string> VH_sig_procs = {"WH_htt","WH_sm_htt", "WH_ps_htt",
            "WH_mm_htt", "ZH_htt","ZH_sm_htt", "ZH_ps_htt",
            "ZH_mm_htt"};        
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
                                             "CMS_eff_trigger_em_13TeV", "lnN", SystMap<>::init(1.02));


        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb,"CMS_eff_t_trg_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb,"CMS_eff_t_trg_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb,"CMS_eff_t_trg_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb,"CMS_eff_t_trg_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb,"CMS_eff_t_trg_MVADM11_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb,"CMS_eff_Xtrigger_mt_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb,"CMS_eff_Xtrigger_mt_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb,"CMS_eff_Xtrigger_mt_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb,"CMS_eff_Xtrigger_mt_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb,"CMS_eff_Xtrigger_mt_MVADM11_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb,"CMS_eff_Xtrigger_et_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb,"CMS_eff_Xtrigger_et_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb,"CMS_eff_Xtrigger_et_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb,"CMS_eff_Xtrigger_et_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).channel({"et_2017","et_2018"}).AddSyst(cb,"CMS_eff_Xtrigger_et_MVADM11_13TeV", "shape", SystMap<>::init(1.00));
        
        //##############################################################################
        //  Electron, muon and tau Id  efficiencies
        //##############################################################################

        cb.cp().AddSyst(cb, "CMS_eff_m", "lnN", SystMap<channel, process>::init
                        ({"mt","mt_2016","mt_2017","mt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, JoinStr({sig_procs, all_mc_bkgs,embed}),  1.01));

        // embedded selection efficiency
        cb.cp().AddSyst(cb, "CMS_eff_m_embedsel", "lnN", SystMap<channel, process>::init
                        ({"em","em_2016","em_2017","em_2018","et","et_2016","et_2017","et_2018","tt","tt_2016","tt_2017","tt_2018","mt","mt_2016","mt_2017","mt_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, embed,  1.04));

        cb.cp().AddSyst(cb, "CMS_eff_e", "lnN", SystMap<channel, process>::init
                        ({"et","et_2016","et_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}, JoinStr({sig_procs, all_mc_bkgs,embed}),       1.02));



        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,4,40}).AddSyst(cb,"CMS_eff_t_pTlow_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,3,30}).AddSyst(cb,"CMS_eff_t_pTlow_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,6,60}).AddSyst(cb,"CMS_eff_t_pTlow_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,5,50}).AddSyst(cb,"CMS_eff_t_pTlow_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).AddSyst(cb,"CMS_eff_t_pTlow_MVADM11_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,4,40}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,3,30}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,6,60}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2,5,50}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}, false).bin_id({1,2}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM11_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM0_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM1_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM2_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM10_13TeV", "shape", SystMap<>::init(1.00));
//        cb.cp().process(JoinStr({sig_procs,all_mc_bkgs,embed})).process({"ZL","EWKZ"},false).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,"CMS_eff_t_pThigh_MVADM11_13TeV", "shape", SystMap<>::init(1.00));


        // the uncorrelated part by channel is due to anti-electron and anti-muon disctiminant, it is the same for MC and embedding so decorrelated in the Morphing. From tau POG recomendation this uncertainty is 3%      

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).process({"ZL","EWKZ"},false).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_mt_13TeV", "lnN", SystMap<>::init(1.03));

        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed})).process({"ZL","EWKZ"},false).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb,
                                             "CMS_eff_t_et_13TeV", "lnN", SystMap<>::init(1.03));

        // 2 real taus
        cb.cp().process(JoinStr({sig_procs, {"ZTT","VVT","TTT","EWKZ"}, embed})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.06));

        // 1 real tau
        cb.cp().process(JoinStr({{"TTJ","ZJ","VVJ","W","Wfakes","TTfakes"}})).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                             "CMS_eff_t_tt_13TeV", "lnN", SystMap<>::init(1.03));       

        //##############################################################################
        //  IP significance correction uncertainty 
        //##############################################################################

        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.000));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.001));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.001));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.000));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"TTT"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.003));
        cb.cp().process({"VVT"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.009));
        cb.cp().process({"Wfakes"}).channel({"tt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.004));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.000));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.002));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.002));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.001));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.000));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.000));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,0.999));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.000));
        cb.cp().process({"TTT"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.980,1.000));
        cb.cp().process({"Wfakes"}).channel({"tt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.001));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.985,1.012));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.005));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.006));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.005));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.005));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,0.999));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.991,1.000));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,0.999));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.989,1.001));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.003));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.005));
        cb.cp().process({"TTT"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.033));
        cb.cp().process({"Wfakes"}).channel({"tt_2016"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.984,0.993));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,1.000));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.990,1.003));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.004));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.002));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.007));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.006));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.006));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.004));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.001));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.001));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.000));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,0.999));
        cb.cp().process({"Wfakes"}).channel({"tt_2016"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.001));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.004));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.007));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.007));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.007));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.005));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.004));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.004));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.001));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.003));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.001));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.002));
        cb.cp().process({"TTT"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.946,1.000));
        cb.cp().process({"ZL"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.949,1.182));
        cb.cp().process({"Wfakes"}).channel({"tt_2016"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.981,1.017));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.988,1.007));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.987,1.013));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.988,1.011));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.005));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.006));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.006));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.997,0.996));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.000));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.989,0.991));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.011));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.006));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(0.996,1.010));
        cb.cp().process({"TTT"}).channel({"tt_2016"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.110,1.000));
        
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.999));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.988));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.989));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.009,0.988));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.988));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.009,0.987));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.988));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.009,0.986));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.988));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.007,0.986));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.011,0.988));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.987));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.989));
        cb.cp().process({"TTT"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.988));
        cb.cp().process({"VVT"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.987));
        cb.cp().process({"ZL"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.996));
        cb.cp().process({"Wfakes"}).channel({"tt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.999));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.989));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.009,0.988));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.989));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.009,0.986));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.988));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.987));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.007,0.986));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.984));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.986));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.987));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.985));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.007,0.985));
        cb.cp().process({"VVT"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.008));
        cb.cp().process({"Wfakes"}).channel({"tt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.998));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.998));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.958));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.036,0.958));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.960));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.960));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.960));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.030,0.960));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.024,0.962));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.023,0.961));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.025,0.962));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.968));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.965));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.968));
        cb.cp().process({"Wfakes"}).channel({"tt_2017"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.005));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.990));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.961));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.962));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.961));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.964));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.966));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.964));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.965));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.965));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.030,0.965));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.024,0.957));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.023,0.956));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.022,0.959));
        cb.cp().process({"VVT"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.104,1.000));
        cb.cp().process({"Wfakes"}).channel({"tt_2017"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.016,1.001));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.998));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.038,0.958));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.039,0.959));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.039,0.958));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.964));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.965));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.964));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.959));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.962));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.961));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.036,0.961));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.038,0.963));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.038,0.962));
        cb.cp().process({"ZL"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.076,0.977));
        cb.cp().process({"Wfakes"}).channel({"tt_2017"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.016,1.001));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.015,1.007));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.072,0.923));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.075,0.927));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.075,0.923));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.059,0.931));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.058,0.936));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.058,0.929));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.068,0.922));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.062,0.925));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.065,0.930));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.086,0.924));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.083,0.922));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.083,0.919));
        cb.cp().process({"ZL"}).channel({"tt_2017"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.15,0.96));        
               
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.999));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.999));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.002));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"TTT"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.006));
        cb.cp().process({"VVT"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.022,1.000));
        cb.cp().process({"ZL"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.937,1.009));
        cb.cp().process({"Wfakes"}).channel({"tt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.001));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.001));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.999));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.004));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.006));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.004));
        cb.cp().process({"TTT"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.999));
        cb.cp().process({"VVT"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.004));
        cb.cp().process({"Wfakes"}).channel({"tt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.995));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,1.000));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,1.000));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,1.000));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.001));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.001));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.994));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.998));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.997));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.994));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.987));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.992));
        cb.cp().process({"VVT"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.082));
        cb.cp().process({"Wfakes"}).channel({"tt_2018"}).bin_id({10}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.972,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.996,0.999));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,0.996));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,0.996));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,0.997));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.997));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.996));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.997));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.994));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.995));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.996));
        cb.cp().process({"TTT"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.928,1.092));
        cb.cp().process({"VVT"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"Wfakes"}).channel({"tt_2018"}).bin_id({9}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.002));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.993));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.999));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.000));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.996));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.995));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.996));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.997));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,0.998));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,0.997));
        cb.cp().process({"TTT"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.050));
        cb.cp().process({"VVT"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.040));
        cb.cp().process({"ZL"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.017));
        cb.cp().process({"Wfakes"}).channel({"tt_2018"}).bin_id({7}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.008));
        cb.cp().process({"EmbedZTT"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.005));
        cb.cp().process({"ggH_sm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.000));
        cb.cp().process({"ggH_ps_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,1.004));
        cb.cp().process({"ggH_mm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"qqH_sm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.007,1.000));
        cb.cp().process({"qqH_ps_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,1.001));
        cb.cp().process({"qqH_mm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.005,1.000));
        cb.cp().process({"WH_sm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.000));
        cb.cp().process({"WH_ps_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.999));
        cb.cp().process({"WH_mm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.000));
        cb.cp().process({"ZH_sm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.011,1.003));
        cb.cp().process({"ZH_ps_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.015,0.998));
        cb.cp().process({"ZH_mm_htt125"}).channel({"tt_2018"}).bin_id({8}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.012,0.999));
       
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.004));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.005));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));

        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.941,1.046));
        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.964,1.041));
        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.948,1.032));
        cb.cp().process({"TTT"}).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.951,1.032)); 

        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.922,1.077));
        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.949,1.132));
        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.932,1.078));
        cb.cp().process({"VVT"}).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.945,1.053));

        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.002));
        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.955,1.037));
        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.965,1.032));
        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.903,1.104));
        cb.cp().process({"ZL"}).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.967,1.030));

        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.006));
        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.003));
        cb.cp().process(ggH_sig_procs).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.003));

        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.004));
        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.003));
        cb.cp().process(qqH_sig_procs).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.003));

        // VH set same as qqH for now
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.000,1.000));
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.997,1.004));
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.003));
        cb.cp().process(VH_sig_procs).channel({"mt_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.998,1.003));

        //////////////////////////////////////////////////////////
        
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.997));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.993));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.997));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.997));

        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.997));
        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.996));
        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.045,0.963));
        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.074,0.935));
        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.026,0.972));
        cb.cp().process({"TTT"}).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.051,0.943)); 

        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.995));
        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.996));
        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.062,0.878));
        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.040,0.968));
        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.952));
        cb.cp().process({"VVT"}).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.018,0.985));

        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.981));
        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.988));
        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.041,0.959));
        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.984));
        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.008,0.990));
        cb.cp().process({"ZL"}).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.007,0.992));

        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.995));
        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.995));
        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.965));
        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.064,0.930));
        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.029,0.968));
        cb.cp().process(ggH_sig_procs).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.968));

        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.995));
        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.994));
        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.968));
        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({4,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.055,0.939));
        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.027,0.968));
        cb.cp().process(qqH_sig_procs).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.968));

        // use qqH numbers also for VH

        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.995));
        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.994));
        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.968));
        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.055,0.939));
        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({5,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.027,0.968));
        cb.cp().process(VH_sig_procs).channel({"mt_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.968));

        //////////////////////////////////////////////////

        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.995));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.011,0.988));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.996));
        cb.cp().process({"EmbedZTT"}).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.994));

        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.049,0.944));
        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.068,0.939));
        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.048,0.960));
        cb.cp().process({"TTT"}).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.043,0.972));

        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.998));
        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998));
        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.025,0.932));
        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.056,0.950));
        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.065,0.949));
        cb.cp().process({"VVT"}).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.933));

        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.996));
        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.026,0.989));
        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.040,0.953));
        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.076,0.927));
        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.969));
        cb.cp().process({"ZL"}).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.979));


        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.997));
        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.997));
        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.020,0.979));
        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.040,0.958));
        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({5,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.018,0.979));
        cb.cp().process(ggH_sig_procs).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.020,0.978));

        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.997));
        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.996));
        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.979));
        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.959));
        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.981));
        cb.cp().process(qqH_sig_procs).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.979));

        // using qqH values also for VH

        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.997));
        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.996));
        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.979));
        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.034,0.959));
        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.981));
        cb.cp().process(VH_sig_procs).channel({"mt_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.979));


        //for ZTT uncertainty is replaced with these shape uncertainties
        cb.cp().channel({"mt_2016"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2016_13TeV", "shape", SystMap<>::init(0.25));
        cb.cp().channel({"mt_2017"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2017_13TeV", "shape", SystMap<>::init(0.25));
        cb.cp().channel({"mt_2018"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2018_13TeV", "shape", SystMap<>::init(0.25));

	////TODO, fix numbers for etau channel

        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.008));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.973,1.024));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.008));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.970,1.029));
        cb.cp().process({"EmbedZTT"}).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_Embed_2018_13TeV", "lnN", SystMapAsymm<>::init(0.965,1.020));

        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.975,1.014));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.959,1.085));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init( 0.992,1.006));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.947,1.049));
        cb.cp().process({"TTT"}).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.929,1.040)); 

        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.987,1.018));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.976,1.019));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.994,1.004));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.969,1.023));
        cb.cp().process({"VVT"}).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.935,1.051));

        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.005));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.008));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.989,1.017));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.991,1.010));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.008));
        cb.cp().process({"ZL"}).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.984,1.034));

        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.997));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.991));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.995));
        cb.cp().process(ggH_sig_procs).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.997));

        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.994));
        cb.cp().process(qqH_sig_procs).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.998));

        // VH set same as qqH for now
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(0.999,1.001));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.994));
        cb.cp().process(VH_sig_procs).channel({"et_2018"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2018_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.998));

        //////////////////////////////////////////////////////////
        
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(0.990,1.008));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.992));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(1.005,0.994));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(0.986,1.011));
        cb.cp().process({"EmbedZTT"}).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_Embed_2017_13TeV", "lnN", SystMapAsymm<>::init(0.986,1.010));

        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.999));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.015,0.985));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.010,0.977));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.006,0.987));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.018,0.982));
        cb.cp().process({"TTT"}).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.969)); 

        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.998));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.024,0.985));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.045,0.942));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.982));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.021,0.993));
        cb.cp().process({"VVT"}).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.024,0.979));

        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.998));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.997));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.041,0.973));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.024,0.971));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.028,0.969));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.048,0.949));
        cb.cp().process({"ZL"}).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.962));

        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.997));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.036,0.957));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.036,0.956));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.021,0.974));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.958));
        cb.cp().process(ggH_sig_procs).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.959));

        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.959));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.960));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.019,0.977));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.958));
        cb.cp().process(qqH_sig_procs).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.958));

        // use qqH numbers also for VH
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.959));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.960));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.019,0.977));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.033,0.958));
        cb.cp().process(VH_sig_procs).channel({"et_2017"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2017_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.958));

        //////////////////////////////////////////////////

        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.993,1.009));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.995,1.003));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(1.003,0.996));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.992,1.006));
        cb.cp().process({"EmbedZTT"}).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_Embed_2016_13TeV", "lnN", SystMapAsymm<>::init(0.991,1.013));

        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.000,0.999));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.021,0.973));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.985));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.991));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.035,0.980));
        cb.cp().process({"TTT"}).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.032,0.973));

        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.021,0.980));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.031,0.960));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.980));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.013,0.990));
        cb.cp().process({"VVT"}).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.004,0.988));

        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.019,0.981));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.976));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.970));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.025,0.979));
        cb.cp().process({"ZL"}).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.022,0.973));

        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.019,0.979));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.021,0.977));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.014,0.982));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({5,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.018,0.982));
        cb.cp().process(ggH_sig_procs).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.020,0.981));

        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.981));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.980));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.011,0.982));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.982));
        cb.cp().process(qqH_sig_procs).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.982));

        // using qqH values also for VH
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({1}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.001,0.999));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({2}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.002,0.998));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({3,30}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.981));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.980));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({4,40}).AddSyst(cb,"CMS_IP_significance_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.011,0.982));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({5,50}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.016,0.982));
        cb.cp().process(VH_sig_procs).channel({"et_2016"}).bin_id({6,60}).AddSyst(cb,"CMS_IP_significance_e_MC_2016_13TeV", "lnN", SystMapAsymm<>::init(1.017,0.982));

        //for ZTT uncertainty is replaced with these shape uncertainties
        cb.cp().channel({"et_2016"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2016_13TeV", "shape", SystMap<>::init(0.25));
        cb.cp().channel({"et_2017"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2017_13TeV", "shape", SystMap<>::init(0.25));
        cb.cp().channel({"et_2018"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2018_13TeV", "shape", SystMap<>::init(0.25));

        //cb.cp().channel({"et_2016"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2016_13TeV", "shapeU", SystMap<>::init(1.0));
        //cb.cp().channel({"et_2017"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2017_13TeV", "shapeU", SystMap<>::init(1.0));
        //cb.cp().channel({"et_2018"}).process({"ZTT"}).AddSyst(cb, "CMS_IPsignifCalib_2018_13TeV", "shapeU", SystMap<>::init(1.0));

        //##############################################################################
        //  b tag and mistag rate  efficiencies 
        //##############################################################################

        cb.cp().channel({"tt","tt_2016","tt_2017","tt_2018"}, false).process({"TTT","VVT",}).AddSyst(cb,
                                             "CMS_eff_b_13TeV", "shape", SystMap<>::init(1.00)); 
        // this is converted to lnN in morphing

        //##############################################################################
        //  Electron, muon and tau energy Scale
        //##############################################################################

        // Add back later!       
 
         cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed, /*{"jetFakes", "QCD"}*/})).channel({"et","et_2016","et_2017","et_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb, 
                                              "CMS_scale_e_13TeV", "shape", SystMap<>::init(1.00)); 
 
       // // Decay Mode based TES Settings
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,4,40}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,3,4,6,30,40,60}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,5,50}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,5}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
        //                                        "CMS_scale_t_3prong1pizero_13TeV", "shape", SystMap<>::init(1.00)); // drop this uncertainty as it is always very small compared to the 3-prong uncertainty

        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,7,8,9,10}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,3,4,5,7,8,9,10,11}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,                                        
                                                "CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,5,6,9,11}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
                                                "CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));
        //cb.cp().process(JoinStr({sig_procs, real_tau_mc_bkgs, embed})).bin_id({1,2,5,6,9,11}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb,
        //                                        "CMS_scale_t_3prong1pizero_13TeV", "shape", SystMap<>::init(1.00)); // drop this uncertainty as it is always very small compared to the 3-prong uncertainty
        

        // Muon 
        // ES Add back after!
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs, embed, /*{"jetFakes", "QCD"}*/})).channel({"mt","mt_2016","mt_2017","mt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                             "CMS_scale_mu_13TeV", "shape", SystMap<>::init(1.00));       
 
        //##############################################################################
        //  Embedded uncertainty on ttbar contamination (and VV contamination)
        //##############################################################################        
        cb.cp().process(embed).AddSyst(cb,"CMS_ttbar_embeded_13TeV", "shape", SystMap<>::init(1.00));
 
        //##############################################################################
        //  jet and met energy Scale
        //##############################################################################
 
        // Add back later!

        // MET Systematic shapes - recoil uncertainties for recoil corrected met, unclustered energy uncertainty for samples with no recoil correction, jes uncertainties propogated to met for samples with no recoil correction
         cb.cp().process({"TT","TTJ","TTT","VV","VVJ","VVT"}).AddSyst(cb, 
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

        cb.cp().process({"ZJ","ZL","ZLL"}).AddSyst(cb,
                                        "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.02));        
 
        cb.cp().process({"ZTT"}).AddSyst(cb,
                                        "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.05));        
 
 
        cb.cp().process({"EWKZ"}).AddSyst(cb,
                                        "CMS_htt_ewkzXsec_13TeV", "lnN", SystMap<>::init(1.05));

        if (! ttbar_fit){
        //   ttbar Normalisation - fully correlated
	    cb.cp().process({"TT","TTT","TTJ"}).AddSyst(cb,
					  "CMS_htt_tjXsec_13TeV", "lnN", SystMap<>::init(1.042));
        }

        // W norm, just for em, tt and the mm region where MC norm is from MC
        
        cb.cp().process({"W"}).channel({"tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018","ttbar","ttbar_2016","ttbar_2017","ttbar_2018"}).AddSyst(cb,
                                                       "CMS_htt_wjXsec_13TeV", "lnN", SystMap<>::init(1.04));
        
        if(control_region==0){
          cb.cp().process({"W"}).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb,
                                                       "CMS_htt_wjXsec_13TeV", "lnN", SystMap<>::init(1.04));    
        }
        
        // QCD uncerts for em (add these if we include the em channel in the combination)
        

        //##############################################################################
        //  Fake-Factor Method Uncertainties
        //##############################################################################

        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_qcd_stat_unc1_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_qcd_stat_unc2_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc1_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,5}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,6}).AddSyst(cb, "ff_et_wjets_stat_unc2_njets2_mvadm2", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_qcd_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_l_pt_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_wjets_syst_njets2", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "ff_et_sub_syst", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_sub_syst", "shape", SystMap<>::init(1.00));

        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc1_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_qcd_stat_unc2_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc1_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,30}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5,50}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6,60}).AddSyst(cb, "ff_mt_wjets_stat_unc2_njets2_mvadm2", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_qcd_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_syst_njets0", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_syst_njets1", "shape", SystMap<>::init(1.00));

        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_l_pt_closure_syst_njets2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_wjets_syst_njets2", "shape", SystMap<>::init(1.00));

        //cb.cp().process({"jetFakes"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).AddSyst(cb, "ff_mt_ttbar_syst", "shape", SystMap<>::init(1.00));

        // split these uncertainties by njets
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_qcd_syst_njets0", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_qcd_syst_njets1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_qcd_met_closure_syst_njets2", "shape", SystMap<>::init(1.00));

        // only need these uncertainties once using the dR corrected FFs
        cb.cp().process({"jetFakes"}).channel({"tt_2016"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2016", "lnN", SystMap<>::init(1.04));
        cb.cp().process({"jetFakes"}).channel({"tt_2017"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2017", "lnN", SystMap<>::init(1.07));
        cb.cp().process({"jetFakes"}).channel({"tt_2018"}).AddSyst(cb, "ff_tt_qcd_syst_norm_2018", "lnN", SystMap<>::init(1.06));

        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_sub_syst", "shape", SystMap<>::init(1.00));


        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc1_njets2_mvadm2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets0_mvadm2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets1_mvadm2", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7,8,9,10}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm0_sig_gt", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm0_sig_lt", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3,4,5,7}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm1", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5,6,9,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm10", "shape", SystMap<>::init(1.00));
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm11", "shape", SystMap<>::init(1.00));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4,10,11}).AddSyst(cb, "ff_tt_qcd_stat_unc2_njets2_mvadm2", "shape", SystMap<>::init(1.00));

        // put an additional uncertainty on the Wfake comonent which is estimated from MC
        cb.cp().process({"Wfakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_tt_mc", "lnN", SystMap<>::init(1.3));

        // additional uncertainty covering non-closures in same-sign data as a function of BDT score
        //cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).AddSyst(cb, "ff_SS_closure", "shape", SystMap<>::init(1.00));
        // this uncertainty would make the ff_tt_mvadm2_closure uncertainty below redundant!

        // additional 5% per sub-leading MVA-DM=2 tau for tt channel to cover non-closures - note this uncertainty is 3% when rounging up for both channels
       
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({4}).AddSyst(cb, "ff_tt_mvadm2_closure", "lnN", SystMap<>::init(1.03));
        cb.cp().process({"jetFakes"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({10}).AddSyst(cb, "ff_tt_mvadm2_closure", "lnN", SystMap<>::init(1.03));

        // for a1-a1 PV method add additional 1% uncertainty to cover SV requirement, this number came from 2016 comparrison so check for 2017 and 2018 also!!

        // prefiring
        cb.cp().process(JoinStr({sig_procs, all_mc_bkgs})).channel({"tt","tt_2016","tt_2017","mt","mt_2016","mt_2017","et","et_2016","et_2017"}).AddSyst(cb,
                                             "CMS_PreFire_13TeV", "shape", SystMap<>::init(1.0));


        //##############################################################################
        //  DY LO->NLO reweighting, Between no and twice the correction.
        //##############################################################################
        
        cb.cp().process( {"ZTT","ZJ","ZL","ZLL"}).AddSyst(cb,
                                             "CMS_htt_dyShape_13TeV", "shape", SystMap<>::init(0.1));
        
        
        //##############################################################################
        // Ttbar shape reweighting, Between no and twice the correction
        //##############################################################################

        
        // add back for mt also
        cb.cp().process( {"TTJ","TTT","TT"}).AddSyst(cb,
                                        "CMS_htt_ttbarShape_13TeV", "shape", SystMap<>::init(1.00));
        
        //##############################################################################
        // ZL shape  and electron/muon  to tau fake only in  mt and et channels (updated March 22)
        //##############################################################################

        // Add back later!
 
        cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,4}).AddSyst(cb,
                                                         "CMS_htt_ZLShape_et_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process( {"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).bin_id({1,2,3,6}).AddSyst(cb,
                                                         "CMS_htt_ZLShape_et_1prong1pi_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb,
                                                         "CMS_htt_ZLShape_mt_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process( {"ZL"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3,6,30,60}).AddSyst(cb,
                                                         "CMS_htt_ZLShape_mt_1prong1pi_13TeV", "shape", SystMap<>::init(1.00));
       
        cb.cp().process({"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "CMS_efake_et_MVADM0_13TeV", "shape", SystMap<>::init(1.0));
        cb.cp().process({"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "CMS_efake_et_MVADM1_13TeV", "shape", SystMap<>::init(1.0));
        cb.cp().process({"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "CMS_efake_et_MVADM2_13TeV", "shape", SystMap<>::init(1.0));
        cb.cp().process({"ZL"}).channel({"et","et_2016","et_2017","et_2018"}).AddSyst(cb, "CMS_efake_et_MVADM10_13TeV", "shape", SystMap<>::init(1.0));

        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.16));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.10));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.14));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.14));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4,40}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_MVADM0_13TeV", "lnN", SystMap<>::init(1.2));

        cb.cp().process({"ZL","EWKZ"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({3,30}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.3));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.05));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.11));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.12));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.08));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.12));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM1_13TeV", "lnN", SystMap<>::init(1.09));

        cb.cp().process({"ZL","EWKZ"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({6,60}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.3));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.015));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.002));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.013));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.013));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.006));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM2_13TeV", "lnN", SystMap<>::init(1.001));

        cb.cp().process({"ZL","EWKZ"}).channel({"mt","mt_2016","mt_2017","mt_2018"}).bin_id({5,50}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.4));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.017));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2016"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.002));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.017));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2017"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.002));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({1}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.018));
        cb.cp().process({"ZL","EWKZ"}).channel({"mt_2018"}).bin_id({2}).AddSyst(cb, "CMS_htt_mFakeTau_MVADM10and11_13TeV", "lnN", SystMap<>::init(1.006));

        cb.cp().process({"ZL","EWKZ"}).channel({"tt","tt_2016"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.10));
        cb.cp().process({"ZL","EWKZ"}).channel({"tt_2017","tt_2018"}).AddSyst(cb,
                                                        "CMS_htt_eFakeTau_13TeV", "lnN", SystMap<>::init(1.02));
        
        //##############################################################################
        // Theoretical Uncertainties on signal (shape/acceptance uncerts needs to be added!)
        //##############################################################################
       
        //scale_gg on signal
        cb.cp().process(JoinStr({ggH_sig_procs, qqH_sig_procs})).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_scale_gg_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({ggH_sig_procs, qqH_sig_procs})).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_PS_ISR_ggH_13TeV", "shape", SystMap<>::init(1.00));

        cb.cp().process(JoinStr({ggH_sig_procs, qqH_sig_procs})).process({"ggH_ph_htt"},false).channel({"et","et_2016","et_2017","et_2018","mt","mt_2016","mt_2017","mt_2018","tt","tt_2016","tt_2017","tt_2018","em","em_2016","em_2017","em_2018"}).AddSyst(cb,
                                             "CMS_PS_FSR_ggH_13TeV", "shape", SystMap<>::init(1.00));
 
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
        cb.cp().process({"WH_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt"}).AddSyst(cb,"QCDscale_VH", "lnN", SystMap<>::init(1.007));
        cb.cp().process({"ZH_htt","ZH_sm_htt","ZH_ps_htt","ZH_mm_htt","ZH_sm_htt","ZH_ps_htt","ZH_mm_htt"}).AddSyst(cb,"QCDscale_VH", "lnN", SystMap<>::init(1.038));
        
        cb.cp().process(JoinStr({ggH_sig_procs, {"ggH_hww125"}})).AddSyst(cb,"pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
        cb.cp().process(JoinStr({qqH_sig_procs, {"qqH_hww125"}})).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
        cb.cp().process({"WH_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt","WH_sm_htt","WH_ps_htt","WH_mm_htt"}).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.019));
        cb.cp().process({"ZH_htt""ZH_sm_htt","ZH_ps_htt","ZH_mm_htt","ZH_sm_htt","ZH_ps_htt","ZH_mm_htt"}).AddSyst(cb,"pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.016));

        
        
    }
}
