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

void Remove13TeVFromNames (ch::CombineHarvester& cb) {
  auto cb_syst = cb.cp();
  cb.cp().ForEachSyst([&](ch::Systematic *syst) {
    std::string old_name = syst->name();
    if (old_name.find("lumi") == std::string::npos) {
      std::string new_name = old_name;
      boost::replace_all(new_name,"_13TeV","");
      syst->set_name(new_name);
    }  
  }
  );
}

void DecorrelateSystSeperateYears (ch::CombineHarvester& cb, string name, std::vector<double> correlations, std::vector<string> chans_2016, std::vector<string> chans_2017, std::vector<string> chans_2018) {
  auto cb_syst = cb.cp().syst_name({name});
  double val2016 = sqrt(1. - correlations[0]);
  double val2017 = sqrt(1. - correlations[1]);
  double val2018 = sqrt(1. - correlations[2]);
  // clone 2016 systs
  if(correlations[0] < 1.) {
    std::cout << "test1" << std::endl;
    ch::CloneSysts(cb.cp().channel(chans_2016).syst_name({name}), cb, [&](ch::Systematic *s) {
        s->set_name(s->name()+"_2016");
        if (s->type().find("shape") != std::string::npos) {
          s->set_scale(s->scale() * val2016);
        }
        if (s->type().find("lnN") != std::string::npos) {
          s->set_value_u((s->value_u() - 1.) * val2016 + 1.);
          if (s->asymm()){
            s->set_value_d((s->value_d() - 1.) * val2016 + 1.);
          }
        }
    });

    if(correlations[0]>0.) {
      // re-scale un-correlated part
      double val = sqrt(correlations[0]);
      cb_syst.channel(chans_2017).ForEachSyst([val](ch::Systematic *syst) {
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
      cb_syst.channel(chans_2016).FilterSysts([&](ch::Systematic *s){
          return s->name().find(name) != std::string::npos && s->name().find("_2016") == std::string::npos && s->name().find("_2017") == std::string::npos && s->name().find("_2018") == std::string::npos;
      });
    }

  }
  // clone 2017 systs
  if(correlations[1] < 1.) {
    ch::CloneSysts(cb.cp().channel(chans_2017).syst_name({name}), cb, [&](ch::Systematic *s) {
        s->set_name(s->name()+"_2017");
        if (s->type().find("shape") != std::string::npos) {
          s->set_scale(s->scale() * val2017);
        }
        if (s->type().find("lnN") != std::string::npos) {
          s->set_value_u((s->value_u() - 1.) * val2017 + 1.);
          if (s->asymm()){
            s->set_value_d((s->value_d() - 1.) * val2017 + 1.);
          }
        }
    });

    if(correlations[1]>0.) {
      // re-scale un-correlated part
      double val = sqrt(correlations[1]);
      cb_syst.channel(chans_2017).ForEachSyst([val](ch::Systematic *syst) {
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
      cb_syst.channel(chans_2017).FilterSysts([&](ch::Systematic *s){
          return s->name().find(name) != std::string::npos && s->name().find("_2016") == std::string::npos && s->name().find("_2017") == std::string::npos && s->name().find("_2018") == std::string::npos;
      });
    }
  }
  // clone 2018 systs
  if(correlations[2] < 1.) {
    ch::CloneSysts(cb.cp().channel(chans_2018).syst_name({name}), cb, [&](ch::Systematic *s) {
        s->set_name(s->name()+"_2018");
        if (s->type().find("shape") != std::string::npos) {
          s->set_scale(s->scale() * val2018);
        }
        if (s->type().find("lnN") != std::string::npos) {
          s->set_value_u((s->value_u() - 1.) * val2018 + 1.);
          if (s->asymm()){
            s->set_value_d((s->value_d() - 1.) * val2018 + 1.);
          }
        }
    });

    if(correlations[2]>0.) {
      // re-scale un-correlated part
      double val = sqrt(correlations[2]);
      cb_syst.channel(chans_2018).ForEachSyst([val](ch::Systematic *syst) {
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
      cb_syst.channel(chans_2018).FilterSysts([&](ch::Systematic *s){
          return s->name().find(name) != std::string::npos && s->name().find("_2016") == std::string::npos && s->name().find("_2017") == std::string::npos && s->name().find("_2018") == std::string::npos;
      });
    }
  }

}



int main(int argc, char** argv) {

    string output_folder = "ztt_validation";
    string input_folder_em="IC/ZTT/";
    string input_folder_et="IC/ZTT/";
    string input_folder_mt="DESY_ZTT/";
    string input_folder_tt="IC/ZTT/";
    string scale_sig_procs="";
    string postfix="";
    unsigned no_shape_systs = 0;
    bool auto_rebin = false;
    bool mergeXbbb = false;
    bool do_embed = false;

    string era;
    po::variables_map vm;
    po::options_description config("configuration");
    config.add_options()
    ("input_folder_em", po::value<string>(&input_folder_em)->default_value("IC/ZTT/"))
    ("input_folder_et", po::value<string>(&input_folder_et)->default_value("IC/ZTT/"))
    ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value("DESY_ZTT/"))
    ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value("IC/ZTT/"))
    ("postfix", po::value<string>(&postfix)->default_value(postfix))
    ("output_folder", po::value<string>(&output_folder)->default_value("ztt_validation"))
    ("no_shape_systs", po::value<unsigned>(&no_shape_systs)->default_value(no_shape_systs))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("era", po::value<string>(&era)->default_value("all"))
    ("mergeXbbb", po::value<bool>(&mergeXbbb)->default_value(false))
    ("do_embed", po::value<bool>(&do_embed)->default_value(false));

    po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
    po::notify(vm);
    typedef vector<string> VString;

    VString years;
    if ( era.find("2016") != std::string::npos ) years.push_back("2016");
    if ( era.find("2017") != std::string::npos ) years.push_back("2017");
    if ( era.find("2018") != std::string::npos ) years.push_back("2018");
    if ( era=="all" ) years = {"2016","2017","2018"};

    typedef vector<pair<int, string>> Categories;

    //! [part1]
    // First define the location of the "auxiliaries" directory where we can
    // source the input files containing the datacard shapes
    //    string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/CombineTools/bin/AllROOT_20fb/";
    //
    std::cout << "input folder = " << input_folder_mt << std::endl;
    std::map<string, string> input_dir;
    input_dir["em"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_em+"/";
    input_dir["mt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_mt+"/";
    input_dir["et"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_et+"/";
    input_dir["tt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_tt+"/";
    input_dir["ttbar"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_em+"/";    
    
    
    VString chns = {"mt"};
    
    map<string, VString> bkg_procs;

    bkg_procs["et"] = {"ZL", "TTT", "VVT", "jetFakes"};
    bkg_procs["mt"] = {"ZL", "TTT", "VVT", "jetFakes"};
    bkg_procs["tt"] = {"ZL", "TTT", "VVT", "jetFakes","Wfakes","TTfakes"};


    ch::CombineHarvester cb;

    
    map<string,Categories> cats;
    
    if( era.find("2016") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["tt_2016"] = {
        {1, "tt_2016_zttEmbed"},
        {2, "tt_2016_jetFakes"}
      };
      cats["mt_2016"] = {
        {3, "mt_murho_alphaLtPiOver4_2016"},
        {4, "mt_mupi_alphaLtPiOver4_2016"},
        {5, "mt_mua1_alphaLtPiOver4_2016"},
        {6, "mt_mu0a1_alphaLtPiOver4_2016"},
        {30, "mt_murho_alphaGtPiOver4_2016"},
        {40, "mt_mupi_alphaGtPiOver4_2016"},
        {50, "mt_mua1_alphaGtPiOver4_2016"},
        {60, "mt_mu0a1_alphaGtPiOver4_2016"},
      };
    } 
    if( era.find("2017") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["tt_2017"] = {
        {3, "tt_2017_ztt_Rho_Rho"},
        //{4, "tt_2017_higgs_0A1_Rho_and_0A1_0A1"},
        {5, "tt_2017_ztt_A1_Rho"},
        {6, "tt_2017_ztt_A1_A1"},
        {7, "tt_2017_ztt_Pi_Rho"},
        {8, "tt_2017_ztt_Pi_Pi"},        
        {9, "tt_2017_ztt_Pi_A1"}, 
        //{10, "tt_2017_higgs_Pi_0A1_Mixed"}, 

      };

      cats["mt_2017"] = {
        {3, "mt_murho_alphaLtPiOver4_2017"},
        {4, "mt_mupi_alphaLtPiOver4_2017"},
        {5, "mt_mua1_alphaLtPiOver4_2017"},
        {6, "mt_mu0a1_alphaLtPiOver4_2017"},
        {30, "mt_murho_alphaGtPiOver4_2017"},
        {40, "mt_mupi_alphaGtPiOver4_2017"},
        {50, "mt_mua1_alphaGtPiOver4_2017"},
        {60, "mt_mu0a1_alphaGtPiOver4_2017"},
      };
    }
    if( era.find("2018") != std::string::npos ||  era.find("all") != std::string::npos) {
      cats["tt_2018"] = {
        {3, "tt_2018_ztt_Rho_Rho"},
        //{4, "tt_2018_higgs_0A1_Rho_and_0A1_0A1"},
        {5, "tt_2018_ztt_A1_Rho"},
        {6, "tt_2018_ztt_A1_A1"},
        {7, "tt_2018_ztt_Pi_Rho"},
        {8, "tt_2018_ztt_Pi_Pi"},
        {9, "tt_2018_ztt_Pi_A1"}, 
        //{10, "tt_2018_higgs_Pi_0A1_Mixed"}, 

      };
      cats["mt_2018"] = {
        {3, "mt_murho_alphaLtPiOver4_2018"},
        {4, "mt_mupi_alphaLtPiOver4_2018"},
        {5, "mt_mua1_alphaLtPiOver4_2018"},
        {6, "mt_mu0a1_alphaLtPiOver4_2018"},
        {30, "mt_murho_alphaGtPiOver4_2018"},
        {40, "mt_mupi_alphaGtPiOver4_2018"},
        {50, "mt_mua1_alphaGtPiOver4_2018"},
        {60, "mt_mu0a1_alphaGtPiOver4_2018"},
      };
    }
    
    map<string, VString> sig_procs;
    sig_procs["ZTT"] = {"ZTT"};
    if(do_embed) sig_procs["ZTT"] = {"EmbedZTT"};

 
    vector<string> masses = {"125"};    
    
    using ch::syst::bin_id;
    
    //! [part2]
    for(auto year: years) {
      for (auto chn : chns) {
          cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, cats[chn+"_"+year]);

          cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, bkg_procs[chn], cats[chn+"_"+year], false);

          if(chn == "em" || chn == "et" || chn == "mt" || chn == "tt"){

            cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn+"_"+year}, sig_procs["ZTT"], cats[chn+"_"+year], true);
          }
      }
    } 
    //! [part4]
    
    
    ch::AddSMRun2Systematics(cb, 0, false, false);
   
    // remove all systematics not relevant for cross check (some could be included in the future for completeness) 
     cb.FilterSysts([&](ch::Systematic *s){
        return s->name().find("eff_b") != std::string::npos || s->name().find("scale_j") != std::string::npos || s->name().find("unclustered") != std::string::npos || s->name().find("ZLShape") != std::string::npos;
      });

      // filter dm=0 syst for now as it appears to be broken - replace with lnN's
      cb.FilterSysts([&](ch::Systematic *s){
        return s->name().find("CMS_eff_t_DM0_13TeV") != std::string::npos;
      });
      using ch::syst::SystMap;
      cb.cp().process({"ZTT","EmbedZTT"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({7,9}).AddSyst(cb,
                                           "CMS_eff_t_DM0_13TeV", "lnN", SystMap<>::init(1.07));
      cb.cp().process({"ZTT","EmbedZTT"}).channel({"tt","tt_2016","tt_2017","tt_2018"}).bin_id({8}).AddSyst(cb,
                                           "CMS_eff_t_DM0_13TeV", "lnN", SystMap<>::init(1.014));

 
    if(no_shape_systs==1){
      cb.FilterSysts([&](ch::Systematic *s){
        return s->type().find("shape") != std::string::npos;
      });
    } else if (no_shape_systs==2){
      // this option will only filter systamtics that required seperate trees to produce. shape systematics made from weights will not be removed
      cb.FilterSysts([&](ch::Systematic *s){
        return s->name().find("scale_t") != std::string::npos || s->name().find("scale_e") != std::string::npos || s->name().find("scale_j") != std::string::npos || s->name().find("_met_") != std::string::npos || s->name().find("ZLShape") != std::string::npos;
      });
    }

    
            
    //! [part7]
    for(auto year: years) {
      for (string chn : chns){
          string channel = chn;
          string extra = "";
          extra = "/"+year+"/";
          cb.cp().channel({chn+"_"+year}).backgrounds().ExtractShapes(
                                                             input_dir[chn] + extra + "htt_"+channel+".inputs-sm-13TeV"+postfix+".root",
                                                             "$BIN/$PROCESS",
                                                             "$BIN/$PROCESS_$SYSTEMATIC");
          if(chn == "em" || chn == "et" || chn == "mt" || chn == "tt"){
            cb.cp().channel({chn+"_"+year}).process(sig_procs["ZTT"]).ExtractShapes(
                                                                    input_dir[chn] + extra + "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                                    "$BIN/$PROCESS",
                                                                    "$BIN/$PROCESS_$SYSTEMATIC");
          }
      }
    }    

    
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
    
    
    
    // can auto-merge the bins with bbb uncertainty > 90% - may be better to merge these bins by hand though!
    auto rebin = ch::AutoRebin()
    .SetBinThreshold(0.)
    .SetBinUncertFraction(0.9)
    .SetRebinMode(1)
    .SetPerformRebin(true)
    .SetVerbosity(1);
    if(auto_rebin) rebin.Rebin(cb, cb);
  
    std::vector<double> binning =  { 0.0000000, 0.62831853, 1.2566371, 1.8849556, 2.5132741, 3.1415927, 3.7699112, 4.3982297, 5.0265482, 5.6548668, 6.2831853 };

    //cb.cp().VariableRebin(binning);
 
  
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

 
    if(mergeXbbb) {
      // if we are mergin bbb's we can't use autoMC stats
      auto bbb = ch::BinByBinFactory()
      .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bbb_bin_$#") // this needs to have "_bbb_bin_" in the pattern for the mergeXbbb option to work
      .SetAddThreshold(0.)
      .SetMergeThreshold(0.4)
      .SetFixNorm(false);
      bbb.MergeBinErrors(cb.cp().backgrounds());
      bbb.AddBinByBin(cb.cp().backgrounds(), cb);


      // if we merge hthe x-axis bins then we need to rename the bbb uncertainties so that they are correlated properly
      // only doing this for the tt channel at the moment, and if we add more channels (no rho-rho) channels for this channel then we might not want to do this for all these categories
      unsigned nxbins=14; // need to hardcode the bin number for the xbins

      cb.cp().backgrounds().channel({"tt","tt_2016","tt_2017","tt_2018"}).ForEachProc([&](ch::Process *proc){
        TH1D *nominal = (TH1D*)proc->ClonedShape().get()->Clone();
        cb.cp().ForEachSyst([&](ch::Systematic *syst) {
          auto old_name = syst->name();
          std::string nonum_name = old_name;
          bool match_proc = (MatchingProcess(*proc,*syst)); 
          if (match_proc && old_name.find("_bbb_bin_") != std::string::npos) {
            int bin_num = -1;
            std::stringstream old_name_ss;
            old_name_ss << old_name;
            string temp;
            int found;
            while (std::getline(old_name_ss, temp, '_')) {
              if (stringstream(temp) >> found) bin_num = found;
            }
            if((bin_num-1) % nxbins==0 ) {
              nonum_name.erase (nonum_name.end()-std::to_string(bin_num).length(), nonum_name.end());
              TH1D *shape_u_new = (TH1D*)syst->ClonedShapeU().get()->Clone();
              TH1D *shape_d_new = (TH1D*)syst->ClonedShapeD().get()->Clone();
              shape_u_new->Add(nominal,-1);
              shape_d_new->Add(nominal,-1);
              std::vector<std::string> names = {};
              for(unsigned i = bin_num+1; i<(unsigned)bin_num+nxbins; ++i) names.push_back(nonum_name+std::to_string(i)); 
              cb.cp().syst_name(names).ForEachSyst([&](ch::Systematic *s) {
                TH1D *shape_u_temp = (TH1D*)s->ClonedShapeU().get()->Clone();
                TH1D *shape_d_temp = (TH1D*)s->ClonedShapeD().get()->Clone();
                shape_u_temp->Add(nominal,-1);
                shape_d_temp->Add(nominal,-1);
                shape_u_new->Add(shape_u_temp);
                shape_d_new->Add(shape_d_temp);
              });
              shape_u_new->Add(nominal);
              shape_d_new->Add(nominal);
              syst->set_shapes(std::unique_ptr<TH1>(static_cast<TH1*>(shape_u_new)),std::unique_ptr<TH1>(static_cast<TH1*>(shape_d_new)),nullptr);
              syst->set_value_u((syst->value_u()-1.)*nxbins + 1.);
              syst->set_value_d((syst->value_d()-1.)*nxbins + 1.); 
              for (auto n : names) {
                cb.FilterSysts([&](ch::Systematic *s){
                  return s->name() == n;
                });
              }
            }  
          }
        });
      });


      // add bbb uncertainties for the signal but as we use reweighted histograms for sm, ps and mm these should be correlated. will need to do something for WH and ZH when we have the samples
      auto bbb_ggh = ch::BinByBinFactory()
      .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_ggH_bin_$#")
      .SetAddThreshold(0.0)
      .SetMergeThreshold(0.0)
      .SetFixNorm(false);
      bbb_ggh.AddBinByBin(cb.cp().signals().process(sig_procs["ggH"]),cb);

      auto bbb_qqh = ch::BinByBinFactory()
      .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_qqH_bin_$#")
      .SetAddThreshold(0.0)
      .SetMergeThreshold(0.0)
      .SetFixNorm(false);
      bbb_qqh.AddBinByBin(cb.cp().signals().process(sig_procs["qqH"]),cb);

    }

    // this part of the code should be used to handle the propper correlations between MC and embedded uncertainties - so no need to try and implement any different treatments in HttSystematics_SMRun2 
  
    // partially decorrelate the energy scale uncertainties
    DecorrelateMCAndEMB(cb,"CMS_scale_e_13TeV","CMS_scale_embedded_e_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_1prong_13TeV","CMS_scale_embedded_t_1prong_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_1prong1pizero_13TeV","CMS_scale_embedded_t_1prong1pizero_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_scale_t_3prong_13TeV","CMS_scale_embedded_t_3prong_13TeV",0.5);
    // partially decorrelate the ID uncertainties uncertainties
    DecorrelateMCAndEMB(cb,"CMS_eff_m","CMS_eff_embedded_m",0.5);
    DecorrelateMCAndEMB(cb,"CMS_eff_e","CMS_eff_embedded_e",0.5);
    DecorrelateMCAndEMB(cb,"CMS_eff_t_mt_13TeV","CMS_eff_embedded_t_mt_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_eff_t_et_13TeV","CMS_eff_embedded_t_et_13TeV",0.5);
    DecorrelateMCAndEMB(cb,"CMS_eff_t_tt_13TeV","CMS_eff_embedded_t_tt_13TeV",0.5);
  
    // fully decorrelate lepton+tau trigger uncertainties for embedded and MC
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_mt_DM0_13TeV","CMS_eff_embedded_Xtrigger_mt_DM0_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_mt_DM1_13TeV","CMS_eff_embedded_Xtrigger_mt_DM1_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_mt_DM10_13TeV","CMS_eff_embedded_Xtrigger_mt_DM10_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_mt_DM11_13TeV","CMS_eff_embedded_Xtrigger_mt_DM11_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_et_DM0_13TeV","CMS_eff_embedded_Xtrigger_et_DM0_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_et_DM1_13TeV","CMS_eff_embedded_Xtrigger_et_DM1_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_et_DM10_13TeV","CMS_eff_embedded_Xtrigger_et_DM10_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_Xtrigger_et_DM11_13TeV","CMS_eff_embedded_Xtrigger_et_DM11_13TeV");
  
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_t_trg_DM0_13TeV","CMS_eff_embedded_t_trg_DM0_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_t_trg_DM1_13TeV","CMS_eff_embedded_t_trg_DM1_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_t_trg_DM10_13TeV","CMS_eff_embedded_t_trg_DM10_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_t_trg_DM11_13TeV","CMS_eff_embedded_t_trg_DM11_13TeV");
  
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_trigger_mt_13TeV","CMS_eff_embedded_trigger_mt_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_trigger_mt_13TeV","CMS_eff_embedded_trigger_et_13TeV");
    cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_eff_trigger_mt_13TeV","CMS_eff_embedded_trigger_em_13TeV");


    // de-correlate systematics for 2016 and 2017, ADD 2018 
    if((era.find("2016") != std::string::npos && era.find("2017") != std::string::npos && era.find("2018") != std::string::npos) ||  era.find("all") != std::string::npos || true){
      std::cout << "Partially Decorrelating systematics for 2016/2017/2018" << std::endl;
      Json::Value js;
      string json_file = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/scripts/correlations.json";
      js = ch::ExtractJsonFromFile(json_file);
      std::vector<std::string> keys = js.getMemberNames();
      for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it){
        string name = *it;
        double value = js[*it].asDouble();
        std::vector<string> chans_2016 = {"em","em_2016","et","et_2016","mt","mt_2016","tt","tt_2016","ttbar","ttbar_2016"};
        std::vector<string> chans_2017 = {"em_2017","et_2017","mt_2017","tt_2017","ttbar_2017"};
        std::vector<string> chans_2018 = {"em_2018","et_2018","mt_2018","tt_2018","ttbar_2018"};
        DecorrelateSyst (cb, name, value, chans_2016, chans_2017, chans_2018);
      }

      // now take care of cases where correlations are different for the 3 years
      string json_file_byyear = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCP2016/scripts/correlations_byyear.json";
      js = ch::ExtractJsonFromFile(json_file_byyear);
      keys = js.getMemberNames();
      for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it){
        string name = *it;
        std::vector<double> values = {js[*it]["2016"].asDouble(),js[*it]["2017"].asDouble(), js[*it]["2018"].asDouble()};
        std::vector<string> chans_2016 = {"em","em_2016","et","et_2016","mt","mt_2016","tt","tt_2016","ttbar","ttbar_2016"};
        std::vector<string> chans_2017 = {"em_2017","et_2017","mt_2017","tt_2017","ttbar_2017"};
        std::vector<string> chans_2018 = {"em_2018","et_2018","mt_2018","tt_2018","ttbar_2018"};
        DecorrelateSystSeperateYears (cb, name, values, chans_2016, chans_2017, chans_2018);
      }

    }

     ch::SetStandardBinNames(cb);
	//! [part8]

	
     // add autoMCStats options
     if(!mergeXbbb) cb.AddDatacardLineAtEnd("* autoMCStats 10 1");
     // add lumi_scale for projection scans
     cb.AddDatacardLineAtEnd("lumi_scale rateParam * *  1. [0,4]");
     cb.AddDatacardLineAtEnd("nuisance edit freeze lumi_scale");

     //! [part9]
     // First we generate a set of bin names:
     
     
     //Write out datacards. Naming convention important for rest of workflow. We
     //make one directory per chn-cat, one per chn and cmb. In this code we only
     //store the individual datacards for each directory to be combined later, but
     //note that it's also possible to write out the full combined card with CH
     string output_prefix = "output/";
     if(output_folder.compare(0,1,"/") == 0) output_prefix="";
     ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$MASS/$BIN.txt",
         	    output_prefix + output_folder + "/$TAG/common/htt_input.root");
     
     
     writer.WriteCards("cmb", cb);
     
     writer.WriteCards("htt_2016", cb.cp().channel({"em_2016","et_2016","mt_2016","tt_2016","ttbar_2016"}));
     writer.WriteCards("htt_2017", cb.cp().channel({"em_2017","et_2017","mt_2017","tt_2017","ttbar_2017"})); 
     writer.WriteCards("htt_2018", cb.cp().channel({"em_2018","et_2018","mt_2018","tt_2018","ttbar_2018"}));
     
     for (auto chn : cb.channel_set()) {
     
       // per-channel
       writer.WriteCards(chn, cb.cp().channel({chn}));
       // And per-channel-category
     }
     //writer.WriteCards("htt_tt_1_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1}));
     //writer.WriteCards("htt_tt_2_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({2}));
     //writer.WriteCards("htt_tt_3_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,3}));
     //writer.WriteCards("htt_tt_4_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,4}));
     //writer.WriteCards("htt_tt_5_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,5}));
     //writer.WriteCards("htt_tt_6_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,6}));
     //writer.WriteCards("htt_tt_7_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,7}));
     //writer.WriteCards("htt_tt_8_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,8}));
     //writer.WriteCards("htt_tt_9_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,9}));
     //writer.WriteCards("htt_tt_10_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,10}));   
     //writer.WriteCards("htt_tt_11_13TeV", cb.cp().channel({"tt_2016","tt_2017","tt_2018"}).bin_id({1,2,11}));
     ///writer.WriteCards("htt_mt_1_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({1}));
     //writer.WriteCards("htt_mt_2_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({2}));
     writer.WriteCards("htt_mt_3_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({3}));
     writer.WriteCards("htt_mt_4_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({4}));
     writer.WriteCards("htt_mt_5_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({5}));
     writer.WriteCards("htt_mt_6_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({6}));

     writer.WriteCards("htt_mt_30_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({30}));
     writer.WriteCards("htt_mt_40_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({40}));
     writer.WriteCards("htt_mt_50_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({50}));
     writer.WriteCards("htt_mt_60_13TeV", cb.cp().channel({"mt_2016","mt_2017","mt_2018"}).bin_id({60}));
        
    cb.PrintAll();
    cout << " done\n";
    
    
}
