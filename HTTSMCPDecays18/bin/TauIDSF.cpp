//***********************************
//
//
//
//Search the word 'ToBeChecked' in this file before running
//
//
//*********************************
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombineTools/interface/CopyTools.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/HTTSMCPDecays18/interface/HttSystematics_SMRun2.h"
#include "CombineHarvester/CombineTools/interface/JsonTools.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"
#include "TF1.h"
#include "TMatrix.h"

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template <typename T>
void To1Bin(T* proc)
{
    std::unique_ptr<TH1> originalHist = proc->ClonedScaledShape();
    TH1F *hist = new TH1F("hist","hist",1,0,1);
    double err = 0;
    double rate =
    originalHist->IntegralAndError(0, originalHist->GetNbinsX() + 1, err);
    hist->SetDirectory(0);
    hist->SetBinContent(1, rate);
    hist->SetBinError(1, err);
    proc->set_shape(*hist, true);  // True means adjust the process rate to the
    // integral of the hist
}



void ConvertShapesToLnN (ch::CombineHarvester& cb, string name, double min_ks) {
  auto cb_syst = cb.cp().syst_name({name});
  cb_syst.ForEachSyst([&](ch::Systematic *syst) {
    if (syst->type().find("shape") != std::string::npos) {
      if(min_ks<=0) {
        std::cout << "Converting systematic " << syst->name() << " for process " << syst->process() << " in bin " << syst->bin() << " to lnN." <<std::endl;
        syst->set_type("lnN");
        return;
      }
      auto shape_u = syst->ClonedShapeU();
      auto shape_d = syst->ClonedShapeD();

      // set uncertainties of up and down templates to zero so they are not concidered in ks test
      for(unsigned i=0; i<=(unsigned)shape_u->GetNbinsX(); ++i) {
        shape_u->SetBinError(i,0.);
        shape_d->SetBinError(i,0.);
      }

      std::unique_ptr<TH1> nominal;

      double ks_u = 0., ks_d = 0.;

      cb.cp().ForEachProc([&](ch::Process *proc){
        bool match_proc = (MatchingProcess(*proc,*syst));
        if(match_proc) nominal = proc->ClonedScaledShape(); 
      });

      if(shape_u && nominal){
        ks_u = nominal.get()->KolmogorovTest(shape_u.get()); 
      }
      if(shape_d && nominal){
        ks_d = nominal.get()->KolmogorovTest(shape_d.get());
      } 
      if(ks_u > min_ks && ks_d > min_ks){
        std::cout << "Converting systematic " << syst->name() << " for process " << syst->process() << " in bin " << syst->bin() << " to lnN. KS scores (u,d) = (" << ks_u << "," << ks_d << ")" <<std::endl;
        syst->set_type("lnN");
      }
      else {
        std::cout << "Not converting systematic " << syst->name() << " for process " << syst->process() << " in bin " << syst->bin() << " to lnN. KS scores (u,d) = (" << ks_u << "," << ks_d << ")" <<std::endl;
      }
    }
  }); 

}

void DecorrelateMCAndEMB (ch::CombineHarvester& cb, string name, string embed_name, double scale) {
  auto cb_syst = cb.cp().process({"EmbedZTT"}).syst_name({name});
  double val = sqrt(1-scale*scale);
  ch::CloneSysts(cb_syst, cb, [&](ch::Systematic *s) {
      s->set_name(embed_name);
      if (s->type().find("shape") != std::string::npos) {
        s->set_scale(s->scale() * val);
      }
      if (s->type().find("lnN") != std::string::npos) {
        s->set_value_u((s->value_u() - 1.) * val + 1.);
        if (s->asymm()){
          s->set_value_d((s->value_d() - 1.) * val + 1.);
        }
      }
  });
  cb_syst.ForEachSyst([scale](ch::Systematic *syst) {
    if (syst->type().find("shape") != std::string::npos) {
      syst->set_scale(syst->scale() * scale);
    }
    if (syst->type().find("lnN") != std::string::npos) {
      syst->set_value_u((syst->value_u() - 1.) * scale + 1.);
      if (syst->asymm()){
        syst->set_value_d((syst->value_d() - 1.) * scale + 1.);
      }
    }
  });

}


void DecorrelateSyst (ch::CombineHarvester& cb, string name, double correlation, std::vector<string> chans_2016, std::vector<string> chans_2017, std::vector<string> chans_2018) {
  if (correlation >= 1.) return;
  auto cb_syst = cb.cp().syst_name({name});
  double val = sqrt(1. - correlation);
  // clone 2016 systs
  ch::CloneSysts(cb.cp().channel(chans_2016).syst_name({name}), cb, [&](ch::Systematic *s) {
      s->set_name(s->name()+"_2016");
      if (s->type().find("shape") != std::string::npos) {
        s->set_scale(s->scale() * val);
      }
      if (s->type().find("lnN") != std::string::npos) {
        s->set_value_u((s->value_u() - 1.) * val + 1.);
        if (s->asymm()){
          s->set_value_d((s->value_d() - 1.) * val + 1.);
        }
      }
  });
  // clone 2017 systs
  ch::CloneSysts(cb.cp().channel(chans_2017).syst_name({name}), cb, [&](ch::Systematic *s) {
      s->set_name(s->name()+"_2017");
      if (s->type().find("shape") != std::string::npos) {
        s->set_scale(s->scale() * val);
      }
      if (s->type().find("lnN") != std::string::npos) {
        s->set_value_u((s->value_u() - 1.) * val + 1.);
        if (s->asymm()){
          s->set_value_d((s->value_d() - 1.) * val + 1.);
        }
      }
  });
  // clone 2018 systs
  ch::CloneSysts(cb.cp().channel(chans_2018).syst_name({name}), cb, [&](ch::Systematic *s) {
      s->set_name(s->name()+"_2018");
      if (s->type().find("shape") != std::string::npos) {
        s->set_scale(s->scale() * val);
      }
      if (s->type().find("lnN") != std::string::npos) {
        s->set_value_u((s->value_u() - 1.) * val + 1.);
        if (s->asymm()){
          s->set_value_d((s->value_d() - 1.) * val + 1.);
        }
      }
  });

  if(correlation>0.) {
    // re-scale un-correlated part
    double val = sqrt(correlation);
    cb_syst.ForEachSyst([val](ch::Systematic *syst) {
      if (syst->type().find("shape") != std::string::npos) {
        syst->set_scale(syst->scale() * val);
      }
      if (syst->type().find("lnN") != std::string::npos) {
        syst->set_value_u((syst->value_u() - 1.) * val + 1.);
        if (syst->asymm()){
          syst->set_value_d((syst->value_d() - 1.) * val + 1.);
        }
      }
    });
  } else {
    // remove uncorrelated part if systs are 100% un-correlated
    cb.FilterSysts([&](ch::Systematic *s){
        return s->name().find(name) != std::string::npos && s->name().find("_2016") == std::string::npos && s->name().find("_2017") == std::string::npos && s->name().find("_2018") == std::string::npos;
    });

  }
}


int main(int argc, char** argv) {


//***************************** ToBeChecked BeforeRunning**************
// A regular command I use to run (specifing these options is recommended): 
// $TauIDSF --era=2016 --embed=true --output_folder='et_datacards_output/embed/2016/'
//
// MC/embed and year will be added to  input folder automatically, but the proper output folder must be given in the option. It just adds 'output/' in the beginning so one should add all the rest (would be good is one make it nicer!)

    //for tight deeptauVsEle (et channel)
    string output_folder = "et_datacards_output/"; //easier to use output option though
    string input_folder_mt="IC/new_et_datacards/"; 
    string input_folder_zmm="IC/et_datacards/";
    
    //for VVLoose deeptauVsEle (mt&tt channel)
    //string output_folder = "ttAndmt_datacards_output/";//easier to use output option though
    //string input_folder_mt="IC/ttAndmt_datacards/";
    //string input_folder_zmm="IC/ttAndmt_datacards/"; 
//***********************************************************************
    
    string scale_sig_procs="";
    string postfix="-2D";
    unsigned no_shape_systs = 0;
    bool auto_rebin = false;
    bool mergeXbbb = false;

    string era;
    bool embed;
    po::variables_map vm;
    po::options_description config("configuration");
    config.add_options()
    //("input_folder_mt", po::value<string>(&input_folder_mt)->default_value("IC/ZTT/"))
    //("input_folder_zmm", po::value<string>(&input_folder_zmm)->default_value("IC/ZTT/"))
    ("postfix", po::value<string>(&postfix)->default_value(postfix))
    ("output_folder", po::value<string>(&output_folder)->default_value("tauIDSF_output"))
    ("no_shape_systs", po::value<unsigned>(&no_shape_systs)->default_value(no_shape_systs))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("era", po::value<string>(&era)->default_value("2018"))
    ("embed", po::value<bool>(&embed)->default_value(false))
    ("mergeXbbb", po::value<bool>(&mergeXbbb)->default_value(false));

    po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
    po::notify(vm);
    typedef vector<string> VString;

    VString years;
    if ( era.find("2016") != std::string::npos ) years.push_back("2016");
    if ( era.find("2017") != std::string::npos ) years.push_back("2017");
    if ( era.find("2018") != std::string::npos ) years.push_back("2018");
    if ( era=="all" ) years = {"2016","2017","2018"};

    typedef vector<pair<int, string>> Categories;

    if (embed){
        input_folder_zmm+="/embed/";
        input_folder_mt+="/embed/";
    } else{
        input_folder_zmm+="/MC/";
        input_folder_mt+="/MC/";
    }


    //! [part1]
    // First define the location of the "auxiliaries" directory where we can
    // source the input files containing the datacard shapes
    //    string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/CombineTools/bin/AllROOT_20fb/";
    std::map<string, string> input_dir;
    input_dir["zmm"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_zmm+"/";
    input_dir["mt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_mt+"/";
    
    
    VString chns = {"zmm","mt"};
    
    map<string, VString> bkg_procs;
    
    if (embed){
        bkg_procs["mt"] = {"ZL","ZJ","TTJ", "TTT","VVJ","VVT","W","QCD"};
        bkg_procs["zmm"] = {"EmbedZL"};
    }
    else{ 
        bkg_procs["mt"] = {"ZL","ZJ","TTJ", "VVJ","W","QCD"};
        bkg_procs["zmm"] = {"ZL", "TT", "VV","W", "ZJ"};
    }


    ch::CombineHarvester cb;
    
    map<string,Categories> cats;
    
    if( era.find("2016") != std::string::npos ||  era.find("all") != std::string::npos) {
     // cats["tt_2016"] = {
       // {2, "tt_2016_jetFakes"}
     // };
     // cats["mt_2016"] = {
      //  {2, "mt_2016_jetFakes"}
     // };
    } 
    if( era.find("2017") != std::string::npos ||  era.find("all") != std::string::npos) {
      //cats["tt_2017"] = {
       // {3, "tt_2017_ztt_Rho_Rho"},
      //};

      //cats["mt_2017"] = {
        //{1, "mt_2018_zttEmbed"},
     // };
    }
   

    if( era.find("2016") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["zmm"] = {
        {1, "zmm_2016_ZMM_inclusive"},
      };
      cats["mt"] = {
        {1, "mt_2016_MVADM0_Pt20to40"},
        {2, "mt_2016_MVADM1_Pt20to40"},
        {3, "mt_2016_MVADM2_Pt20to40"},
        {4, "mt_2016_MVADM10_Pt20to40"},
        {5, "mt_2016_MVADM11_Pt20to40"},
        {6, "mt_2016_MVADM0_PtMoreThan40"},
        {7, "mt_2016_MVADM1_PtMoreThan40"},
        {8, "mt_2016_MVADM2_PtMoreThan40"},
        {9, "mt_2016_MVADM10_PtMoreThan40"},
        {10, "mt_2016_MVADM11_PtMoreThan40"},
        
        {11, "mt_2016_HPSDM0_Pt20to40"},
        {12, "mt_2016_HPSDM1_Pt20to40"},
        {13, "mt_2016_HPSDM10_Pt20to40"},
        {14, "mt_2016_HPSDM11_Pt20to40"},
        {15, "mt_2016_HPSDM0_PtMoreThan40"},
        {16, "mt_2016_HPSDM1_PtMoreThan40"},
        {17, "mt_2016_HPSDM10_PtMoreThan40"},
        {18, "mt_2016_HPSDM11_PtMoreThan40"},
        
        {19, "mt_2016_MVADM1_NoHPS0_Pt20to40"},
        {20, "mt_2016_MVADM2_NoHPS0_Pt20to40"},
        {21, "mt_2016_MVADM1_NoHPS0_PtMoreThan40"},
        {22, "mt_2016_MVADM2_NoHPS0_PtMoreThan40"},
      };
      
    }


    if( era.find("2017") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["zmm"] = {
        {1, "zmm_2017_ZMM_inclusive"},
      };
      cats["mt"] = {
        {1, "mt_2017_MVADM0_Pt20to40"},
        {2, "mt_2017_MVADM1_Pt20to40"},
        {3, "mt_2017_MVADM2_Pt20to40"},
        {4, "mt_2017_MVADM10_Pt20to40"},
        {5, "mt_2017_MVADM11_Pt20to40"},
        {6, "mt_2017_MVADM0_PtMoreThan40"},
        {7, "mt_2017_MVADM1_PtMoreThan40"},
        {8, "mt_2017_MVADM2_PtMoreThan40"},
        {9, "mt_2017_MVADM10_PtMoreThan40"},
        {10, "mt_2017_MVADM11_PtMoreThan40"},
        
        {11, "mt_2017_HPSDM0_Pt20to40"},
        {12, "mt_2017_HPSDM1_Pt20to40"},
        {13, "mt_2017_HPSDM10_Pt20to40"},
        {14, "mt_2017_HPSDM11_Pt20to40"},
        {15, "mt_2017_HPSDM0_PtMoreThan40"},
        {16, "mt_2017_HPSDM1_PtMoreThan40"},
        {17, "mt_2017_HPSDM10_PtMoreThan40"},
        {18, "mt_2017_HPSDM11_PtMoreThan40"},
        
        {19, "mt_2017_MVADM1_NoHPS0_Pt20to40"},
        {20, "mt_2017_MVADM2_NoHPS0_Pt20to40"},
        {21, "mt_2017_MVADM1_NoHPS0_PtMoreThan40"},
        {22, "mt_2017_MVADM2_NoHPS0_PtMoreThan40"},
      };
    }


    if( era.find("2018") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["zmm"] = {
        {1, "zmm_2018_ZMM_inclusive"},
      };
      cats["mt"] = {
        {1, "mt_2018_MVADM0_Pt20to40"},
        {2, "mt_2018_MVADM1_Pt20to40"},
        {3, "mt_2018_MVADM2_Pt20to40"},
        {4, "mt_2018_MVADM10_Pt20to40"},
        {5, "mt_2018_MVADM11_Pt20to40"},
        {6, "mt_2018_MVADM0_PtMoreThan40"},
        {7, "mt_2018_MVADM1_PtMoreThan40"},
        {8, "mt_2018_MVADM2_PtMoreThan40"},
        {9, "mt_2018_MVADM10_PtMoreThan40"},
        {10, "mt_2018_MVADM11_PtMoreThan40"},
        
        {11, "mt_2018_HPSDM0_Pt20to40"},
        {12, "mt_2018_HPSDM1_Pt20to40"},
        {13, "mt_2018_HPSDM10_Pt20to40"},
        {14, "mt_2018_HPSDM11_Pt20to40"},
        {15, "mt_2018_HPSDM0_PtMoreThan40"},
        {16, "mt_2018_HPSDM1_PtMoreThan40"},
        {17, "mt_2018_HPSDM10_PtMoreThan40"},
        {18, "mt_2018_HPSDM11_PtMoreThan40"},
        
        {19, "mt_2018_MVADM1_NoHPS0_Pt20to40"},
        {20, "mt_2018_MVADM2_NoHPS0_Pt20to40"},
        {21, "mt_2018_MVADM1_NoHPS0_PtMoreThan40"},
        {22, "mt_2018_MVADM2_NoHPS0_PtMoreThan40"},
      };
    }


    map<string, VString> sig_procs;
    
    if (embed) sig_procs["mt"] = {"EmbedZTT"};
    else sig_procs["mt"] = {"ZTT","TTT", "VVT"};
    
    vector<string> masses = {"125"};    
    
    using ch::syst::bin_id;
    
    //! [part2]
    for (auto chn : chns) {
        cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn}, cats[chn]);

        cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {chn}, bkg_procs[chn], cats[chn], false);

        if( chn == "mt"){

          cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn}, sig_procs["mt"], cats[chn], true);
        }
    }
    //! [part4]
    
    
    // Add systematics here  
      std::vector<std::string> all_mc = {"ZL","ZJ","ZTT","TTJ","TTT","TT","W","VV","VVT","VVJ"};
      std::vector<std::string> all_mc_noW = {"ZL","ZJ","ZTT","TTJ","TTT","TT","VV","VVT","VVJ"};
      std::vector<std::string> all_mc_T = {"ZTT","TTT","VVT"};
      std::vector<std::string> embed_proc = {"EmbedZL","EmbedZTT"};
      
      using ch::syst::SystMap;
      using ch::syst::process;
      using ch::syst::channel;
      using ch::JoinStr;

      if (embed) cb.cp().process({"ZL","ZJ","ZTT","TTJ","TTT","TT","VV","VVT","VVJ"}).AddSyst(cb, "lumi_13TeV", "lnN", SystMap<>::init(1.025));
      else cb.cp().process({"TTJ","TTT","TT","VV","VVT","VVJ"}).AddSyst(cb, "lumi_13TeV", "lnN", SystMap<>::init(1.025));


      if (embed) cb.cp().process({"EmbedZL","EmbedZTT"}).AddSyst(cb, "rate_Z", "rateParam", SystMap<>::init(1.00));   
      else cb.cp().process({"ZTT","ZL","ZJ"}).AddSyst(cb, "rate_Z", "rateParam", SystMap<>::init(1.00));
           

        //##############################################################################
        //  trigger   
        //##############################################################################

        if(embed) {
          cb.cp().AddSyst(cb, "CMS_eff_trigger_mt_13TeV", "lnN", SystMap<channel, process>::init
                          ({"mt"}, {"TTT","VVT","ZL","ZJ","TTJ","VVJ"},  1.02)
                          );
        } else {

          cb.cp().AddSyst(cb, "CMS_eff_trigger_mt_13TeV", "lnN", SystMap<channel, process>::init
                          ({"mt"}, {"TTT","VVT","TTJ","VVJ"},  1.02)
                          ({"zmm"}, {"TT","VV"},  1.02)
                          );

        }


        //##############################################################################
        //  Electron, muon and tau Id  efficiencies
        //##############################################################################

        if(embed) { 
          cb.cp().AddSyst(cb, "CMS_eff_m", "lnN", SystMap<channel, process>::init
                          ({"mt"}, {embed_proc},  0.99)
                          ({"mt"}, {"TTT","VVT","ZL","ZJ","TTJ","VVJ"},  1.01)
                          );
        } else {

          cb.cp().AddSyst(cb, "CMS_eff_m", "lnN", SystMap<channel, process>::init
                          ({"mt"}, {"ZTT","ZL","ZJ"},  0.99)
                          ({"mt"}, {"TTT","VVT","TTJ","VVJ"},  1.01)
                          ({"zmm"}, {"TT","VV"},  1.02)
                          );

        }

        // embedded selection efficiency
        //cb.cp().process(embed_proc).AddSyst(cb,
        //                                     "CMS_eff_m_embedsel", "lnN", SystMap<>::init(1.04));

        //##############################################################################
        //   muon and tau energy Scale
        //##############################################################################

       // // Decay Mode based TES Settings
        cb.cp().process(JoinStr({all_mc_T,{"EmbedZTT"}})).channel({"mt"}).AddSyst(cb,
                                      "CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({all_mc_T,{"EmbedZTT"}})).channel({"mt"}).AddSyst(cb,
                                      "CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({all_mc_T,{"EmbedZTT"}})).channel({"mt"}).AddSyst(cb,
                                      "CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));
        cb.cp().process(JoinStr({all_mc_T,{"EmbedZTT"}})).channel({"mt"}).AddSyst(cb,
                                      "CMS_scale_t_3prong1pizero_13TeV", "shape", SystMap<>::init(1.00));

        // Muon 
        cb.cp().process(JoinStr({ all_mc,embed_proc})).channel({"mt","zmm"}).AddSyst(cb,
                                             "CMS_scale_mu_13TeV", "shape", SystMap<>::init(1.00));       
 




        //##############################################################################
        //  jet and met energy Scale
        //##############################################################################
 

         cb.cp().process({"ZTT","ZL","ZJ"}).channel({"mt"}).AddSyst(cb,
                                                   "CMS_htt_boson_reso_met_13TeV", "shape", SystMap<>::init(1.00));  
         cb.cp().process({"ZTT","ZL","ZJ"}).channel({"mt"}).AddSyst(cb,
                                                   "CMS_htt_boson_scale_met_13TeV", "shape", SystMap<>::init(1.00));      



        //##############################################################################
        //  Background normalization uncertainties
        //##############################################################################
        
        //   Diboson  Normalisation - fully correlated
        cb.cp().process({"VV","VVT","VVJ"}).AddSyst(cb,
                                        "CMS_htt_vvXsec_13TeV", "lnN", SystMap<>::init(1.05));

        if(embed) {  
          cb.cp().process({"ZTT","ZJ","ZL"}).AddSyst(cb,
                                          "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.02));        
        }

        //   ttbar Normalisation - fully correlated
	    cb.cp().process({"TT","TTT","TTJ"}).AddSyst(cb,
					  "CMS_htt_tjXsec_13TeV", "lnN", SystMap<>::init(1.042));


        //##############################################################################
        //  DY LO->NLO reweighting, Between no and twice the correction.
        //##############################################################################
        
        cb.cp().process( {"ZTT","ZJ","ZL"}).AddSyst(cb,
                                             "CMS_htt_dyShape_13TeV", "shape", SystMap<>::init(0.1));
        
        
        //##############################################################################
        // Ttbar shape reweighting, Between no and twice the correction
        //##############################################################################
        
        cb.cp().process( {"TTJ","TTT","TT"}).AddSyst(cb,
                                        "CMS_htt_ttbarShape_13TeV", "shape", SystMap<>::init(1.00));


        // weighted avarages of recommended tau POG uncertainties provided in bins of eta (update later!)
        cb.cp().process({"ZL"}).channel({"mt"}).AddSyst(cb,
                                                        "CMS_htt_mFakeTau_13TeV", "lnN", SystMap<>::init(1.4));


        cb.cp().process({"VVJ","TTJ"}).channel({"mt"}).AddSyst(cb,
                                                        "CMS_htt_JFakeTau_13TeV", "lnN", SystMap<>::init(1.2));
    
    
        cb.cp().process({"W"}).channel({"mt"}).AddSyst(cb,
                                                        "CMS_htt_WYield_13TeV", "lnN", SystMap<>::init(1.1));
    
        cb.cp().process({"QCD"}).channel({"mt"}).AddSyst(cb,
                                                        "CMS_htt_QCDYield_13TeV", "lnN", SystMap<>::init(1.1));
    
    
        if (embed) cb.cp().process({"TTT","VVT"}).channel({"mt"}).AddSyst(cb, "CMS_htt_TauID_13TeV", "lnN", SystMap<>::init(1.15));


    //! [part7]
    for(auto year: years) {
      for (string chn : chns){
          string channel = chn;
          string extra = "";
          extra = "/"+year+"/";
          cb.cp().channel({chn}).backgrounds().ExtractShapes(
                                                             input_dir[chn] + extra + "htt_"+channel+".inputs-sm-13TeV"+postfix+".root",
                                                             "$BIN/$PROCESS",
                                                             "$BIN/$PROCESS_$SYSTEMATIC");
          if(chn == "mt"){
            cb.cp().channel({chn}).process(sig_procs["mt"]).ExtractShapes(
                                                                    input_dir[chn] + extra + "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                                    "$BIN/$PROCESS",
                                                                    "$BIN/$PROCESS_$SYSTEMATIC");
          }
      }
    }    

//convert shape to lnN
        cb.cp().channel({"zmm"}).syst_type({"shape"}).ForEachSyst([](ch::Systematic *sys) {sys->set_type("lnN"); });

        ConvertShapesToLnN(cb.cp().process({"ZTT","ZL","ZJ"}), "CMS_htt_boson_reso_met_13TeV", 0.);
        ConvertShapesToLnN(cb.cp().process({"ZTT","ZL","ZJ"}), "CMS_htt_boson_scale_met_13TeV", 0.);
    
    //Now delete processes with 0 yield
    cb.FilterProcs([&](ch::Process *p) {
        bool null_yield = !(p->rate() > 0.);
        if (null_yield){
            std::cout << "[Null yield] Removing process with null yield: \n ";
            std::cout << ch::Process::PrintHeader << *p << "\n";
            cb.FilterSysts([&](ch::Systematic *s){
                bool remove_syst = (MatchingProcess(*p,*s));
                return remove_syst;
            });
        }
        return null_yield;
    });   
    
  
    // At this point we can fix the negative bins
    std::cout << "Fixing negative bins\n";
    cb.ForEachProc([](ch::Process *p) {
      if (ch::HasNegativeBins(p->shape())) {
         std::cout << "[Negative bins] Fixing negative bins for " << p->bin()
                   << "," << p->process() << "\n";
        auto newhist = p->ClonedShape();
        ch::ZeroNegativeBins(newhist.get());
        // Set the new shape but do not change the rate, we want the rate to still
        // reflect the total integral of the events
        p->set_shape(std::move(newhist), false);
      }
    });
  
    cb.ForEachSyst([](ch::Systematic *s) {
      if (s->type().find("shape") == std::string::npos) return;            
      if (ch::HasNegativeBins(s->shape_u()) || ch::HasNegativeBins(s->shape_d())) {
         std::cout << "[Negative bins] Fixing negative bins for syst" << s->bin()
               << "," << s->process() << "," << s->name() << "\n";
        auto newhist_u = s->ClonedShapeU();
        auto newhist_d = s->ClonedShapeD();
        ch::ZeroNegativeBins(newhist_u.get());
        ch::ZeroNegativeBins(newhist_d.get());
        // Set the new shape but do not change the rate, we want the rate to still
        // reflect the total integral of the events
        s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
      }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
      if (s->type().find("shape") == std::string::npos) return;
      if(!(s->value_d()<0.001 || s->value_u()<0.001)) return;
      std::cout << "[Negative yield] Fixing negative yield for syst" << s->bin()
               << "," << s->process() << "," << s->name() << "\n";
      if(s->value_u()<0.001){
         s->set_value_u(0.001);
         s->set_type("lnN");
      }
      if(s->value_d()<0.001){
         s->set_value_d(0.001);
         s->set_type("lnN");
      }
  });



    // this part of the code should be used to handle the propper correlations between MC and embedded uncertainties - so no need to try and implement any different treatments in HttSystematics_SMRun2 
  
    // partially decorrelate the energy scale uncertainties
    DecorrelateMCAndEMB(cb,"CMS_scale_t_1prong_13TeV","CMS_scale_embedded_t_1prong_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_1prong1pizero_13TeV","CMS_scale_embedded_t_1prong1pizero_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_3prong_13TeV","CMS_scale_embedded_t_3prong_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_3prong1pizero_13TeV","CMS_scale_embedded_t_3prong1pi0_13TeV",0.5);
    // partially decorrelate the ID uncertainties uncertainties
    DecorrelateMCAndEMB(cb,"CMS_eff_m","CMS_eff_embedded_m",0.5);
  
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_trigger_mt_13TeV","CMS_eff_embedded_trigger_mt_13TeV");


     ch::SetStandardBinNames(cb);
	//! [part8]

	
     // add autoMCStats options
     //cb.AddDatacardLineAtEnd("* autoMCStats 10 1");
     cb.AddDatacardLineAtEnd("* autoMCStats 0 1");
     // add lumi_scale for projection scans
     cb.AddDatacardLineAtEnd("lumi_scale rateParam * *  1. [0,4]");
     cb.AddDatacardLineAtEnd("nuisance edit freeze lumi_scale");

     // define categories 
    cb.cp().channel({"zmm"}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","zmm_cat");}); 
    // copy the line below and add cati for every mt category
    cb.cp().channel({"mt"}).bin_id({1}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat1");});
    cb.cp().channel({"mt"}).bin_id({2}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat2");});
    cb.cp().channel({"mt"}).bin_id({3}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat3");});
    cb.cp().channel({"mt"}).bin_id({4}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat4");});
    cb.cp().channel({"mt"}).bin_id({5}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat5");});
    cb.cp().channel({"mt"}).bin_id({6}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat6");});
    cb.cp().channel({"mt"}).bin_id({7}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat7");});
    cb.cp().channel({"mt"}).bin_id({8}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat8");});
    cb.cp().channel({"mt"}).bin_id({9}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat9");});
    cb.cp().channel({"mt"}).bin_id({10}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat10");});

    cb.cp().channel({"mt"}).bin_id({11}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat11");});
    cb.cp().channel({"mt"}).bin_id({12}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat12");});
    cb.cp().channel({"mt"}).bin_id({13}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat13");});
    cb.cp().channel({"mt"}).bin_id({14}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat14");});
    cb.cp().channel({"mt"}).bin_id({15}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat15");});
    cb.cp().channel({"mt"}).bin_id({16}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat16");});
    cb.cp().channel({"mt"}).bin_id({17}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat17");});
    cb.cp().channel({"mt"}).bin_id({18}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat18");});
    
    cb.cp().channel({"mt"}).bin_id({19}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat19");});
    cb.cp().channel({"mt"}).bin_id({20}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat20");});
    cb.cp().channel({"mt"}).bin_id({21}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat21");});
    cb.cp().channel({"mt"}).bin_id({22}).ForEachObj([&](ch::Object *obj){obj->set_attribute("cat","cat22");});
    
    //! [part9]
     // First we generate a set of bin names:
     
     
     //Write out datacards. Naming convention important for rest of workflow. We
     //make one directory per chn-cat, one per chn and cmb. In this code we only
     //store the individual datacards for each directory to be combined later, but
     //note that it's also possible to write out the full combined card with CH
     string output_prefix = "output/";
     if(output_folder.compare(0,1,"/") == 0) output_prefix="";
     std::cout<<"out_dir="<<output_prefix + output_folder + "/$TAG/$MASS/$BIN.txt"<<std::endl;
     ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$MASS/$BIN.txt",
         	    output_prefix + output_folder + "/$TAG/common/htt_input.root");
      
     writer.WriteCards("cmb", cb);
     writer.WriteCards("htt_mt_1_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat1","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_2_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat2","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_3_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat3","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_4_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat4","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_5_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat5","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_6_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat6","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_7_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat7","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_8_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat8","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_9_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat9","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_10_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat10","zmm_cat"},"cat"));

     writer.WriteCards("htt_mt_11_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat11","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_12_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat12","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_13_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat13","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_14_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat14","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_15_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat15","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_16_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat16","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_17_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat17","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_18_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat18","zmm_cat"},"cat"));
     
     writer.WriteCards("htt_mt_19_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat19","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_20_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat20","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_21_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat21","zmm_cat"},"cat"));
     writer.WriteCards("htt_mt_22_13TeV", cb.cp().channel({"mt","zmm"}).attr({"cat22","zmm_cat"},"cat"));
     
     //for (auto chn : cb.channel_set()) {
     //
     //  // per-channel
     //  writer.WriteCards(chn, cb.cp().channel({chn}));
     //  // And per-channel-category
     //}
     //writer.WriteCards("htt_tt_1_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1}));
     //writer.WriteCards("htt_tt_2_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({2}));
    // writer.WriteCards("htt_tt_3_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3}));
    // //writer.WriteCards("htt_tt_4_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4}));
    // writer.WriteCards("htt_tt_5_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5}));
    // writer.WriteCards("htt_tt_6_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,6}));
    // writer.WriteCards("htt_tt_7_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7}));
    // writer.WriteCards("htt_tt_8_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,8}));
    // writer.WriteCards("htt_tt_9_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,9}));
    // //writer.WriteCards("htt_tt_10_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,10}));   
    // //writer.WriteCards("htt_tt_11_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,11}));
    // //writer.WriteCards("htt_mt_1_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1}));
    // //writer.WriteCards("htt_mt_2_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({2}));
    // writer.WriteCards("htt_mt_3_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1,2,3}));
    // writer.WriteCards("htt_mt_4_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1,2,4}));
    // writer.WriteCards("htt_mt_5_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1,2,5}));
     //writer.WriteCards("htt_mt_6_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1,2,6}));
        
    cb.PrintAll();
    cout << " done\n";
    
    
}
