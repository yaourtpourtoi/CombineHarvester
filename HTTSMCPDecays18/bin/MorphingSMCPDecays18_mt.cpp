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

void DecorrelateSyst (ch::CombineHarvester& cb, string name, double correlation, std::vector<string> chans_2016, std::vector<string> chans_2017) {
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
        return s->name().find(name) != std::string::npos && s->name().find("_2016") == std::string::npos && s->name().find("_2017") == std::string::npos;
    });

  }
}

int main(int argc, char** argv) {

    string output_folder = "sm_run2";
    string input_folder_em="Imperial/CP/";
    string input_folder_et="Imperial/CP/";
    string input_folder_mt="Imperial/CP/";
    string input_folder_tt="Imperial/CP/";
    string input_folder_mm="USCMS/";
    string scale_sig_procs="";
    string postfix="";
    bool ttbar_fit = false;
    unsigned no_shape_systs = 0;
    bool do_embedding = true;
    bool auto_rebin = false;
    bool do_jetfakes = true;
    bool do_mva = false;    
    int do_control_plots = 0;
    bool doDecays = false;
    bool mergeXbbb = false; 
    bool id_cats = false;

    string era;
    po::variables_map vm;
    po::options_description config("configuration");
    config.add_options()
    ("input_folder_em", po::value<string>(&input_folder_em)->default_value("Imperial/CP"))
    ("input_folder_et", po::value<string>(&input_folder_et)->default_value("Imperial/CP"))
    ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value("Imperial/CP"))
    ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value("Imperial/CP"))
    ("input_folder_mm", po::value<string>(&input_folder_mm)->default_value("USCMS"))
    ("postfix", po::value<string>(&postfix)->default_value(postfix))
    ("output_folder", po::value<string>(&output_folder)->default_value("sm_run2"))
    ("no_shape_systs", po::value<unsigned>(&no_shape_systs)->default_value(no_shape_systs))
    ("do_embedding", po::value<bool>(&do_embedding)->default_value(true))
    ("do_jetfakes", po::value<bool>(&do_jetfakes)->default_value(true))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("do_mva", po::value<bool>(&do_mva)->default_value(false))
    ("do_control_plots", po::value<int>(&do_control_plots)->default_value(0))    
    ("era", po::value<string>(&era)->default_value("2016"))
    ("ttbar_fit", po::value<bool>(&ttbar_fit)->default_value(false))
    ("doDecays", po::value<bool>(&doDecays)->default_value(false))
    ("mergeXbbb", po::value<bool>(&mergeXbbb)->default_value(false))
    ("id_cats", po::value<bool>(&id_cats)->default_value(false));

    po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
    po::notify(vm);
    typedef vector<string> VString;

 
    if(do_control_plots>0){
      ttbar_fit = false;
      input_folder_em="/Imperial/control_cards_"+era+"/";
      input_folder_et="/Imperial/control_cards_"+era+"/";
      input_folder_mt="/Imperial/control_cards_"+era+"/";
      input_folder_tt="/Imperial/control_cards_"+era+"/";
    }

 
    VString years;
    if ( era.find("2016") != std::string::npos ) years.push_back("2016");
    if ( era.find("2017") != std::string::npos ) years.push_back("2017");
    if ( era=="all" ) years = {"2016","2017"};
 
    typedef vector<string> VString;
    typedef vector<pair<int, string>> Categories;
    //! [part1]
    // First define the location of the "auxiliaries" directory where we can
    // source the input files containing the datacard shapes
    //    string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/CombineTools/bin/AllROOT_20fb/";
    std::map<string, string> input_dir;
    input_dir["em"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_em+"/";
    input_dir["mt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_mt+"/";
    input_dir["et"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_et+"/";
    input_dir["tt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_tt+"/";
    input_dir["ttbar"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSMCPDecays18/shapes/"+input_folder_em+"/";    
    
    
    VString chns = {"tt"};
    if (ttbar_fit) chns.push_back("ttbar");
    
    map<string, VString> bkg_procs;
    bkg_procs["et"] = {"ZTT", "QCD", "ZL", "ZJ","TTT","TTJ", "VVT", "VVJ", "EWKZ", "W"};
    bkg_procs["mt"] = {"ZTT", "QCD", "ZL", "ZJ","TTT","TTJ", "VVT", "VVJ", "EWKZ", "W"};
    bkg_procs["tt"] = {"ZTT", "W", "QCD", "ZL", "ZJ","TTT","TTJ",  "VVT","VVJ", "EWKZ"};
    bkg_procs["em"] = {"ZTT","W", "QCD", "ZLL", "TT", "VV", "EWKZ"};
    bkg_procs["ttbar"] = {"ZTT", "W", "QCD", "ZLL", "TT", "VV", "EWKZ"};
    
    if(do_embedding){
      bkg_procs["et"] = {"EmbedZTT", "QCD", "ZL", "ZJ","TTT","TTJ", "VVT", "VVJ", "W", "EWKZ"};
      bkg_procs["mt"] = {"EmbedZTT", "QCD", "ZL", "ZJ","TTT","TTJ",  "VVT", "VVJ", "W", "EWKZ"};
      bkg_procs["tt"] = {"EmbedZTT", "W", "QCD", "ZL", "ZJ","TTT","TTJ",  "VVT","VVJ", "EWKZ"};
      bkg_procs["em"] = {"EmbedZTT","W", "QCD", "ZLL", "TT", "VV", "EWKZ"};
      bkg_procs["ttbar"] = {"EmbedZTT", "W", "QCD", "ZLL", "TT", "VV", "EWKZ"};
    }

    if(do_jetfakes){
      bkg_procs["et"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes"};
      bkg_procs["mt"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes"};
      bkg_procs["tt"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes"};

      if(do_embedding){
        bkg_procs["et"] = {"EmbedZTT", "ZL", "TTT", "VVT", "jetFakes", "EWKZ"};
        bkg_procs["mt"] = {"EmbedZTT", "ZL", "TTT", "VVT", "jetFakes", "EWKZ"};
        bkg_procs["tt"] = {"EmbedZTT", "ZL", "TTT", "VVT", "jetFakes", "EWKZ"};
      }

    }


    ch::CombineHarvester cb;
    
    map<string,Categories> cats;
    
    map<string,Categories> cats_cp;

    if (doDecays && !do_mva) {
      cats_cp["tt_2016"] = {
        {1, "tt_0jet_rho"},
        {2, "tt_boosted_rho"},
        //{3, "tt_dijet_rho"}
        {3, "tt_dijet_lowboost_rho"},
        {4, "tt_dijet_boosted_rho"}

        /* {5, "tt_0jet_other"}, */
        /* {6, "tt_boosted_other"}, */
        /* //{3, "tt_dijet_rho"} */
        /* {7, "tt_dijet_lowboost_other"}, */
        /* {8, "tt_dijet_boosted_other"} */

        /*{1, "tt_0jet_rho_idg0p5"},
        {2, "tt_boosted_rho_idg0p5"},
        {3, "tt_dijet_rho_idg0p5"},

        {4, "tt_0jet_rho_!idg0p5"},
        {5, "tt_boosted_rho_!idg0p5"},
        {6, "tt_dijet_rho_!idg0p5"}*/
      };

      cats_cp["mt_2016"] = {
        {1, "mt_0jet_mixed"},
        {2, "mt_boosted_mixed"},
        //{3, "mt_dijet_mixed"}
        {3, "mt_dijet_lowboost_mixed"},
        {4, "mt_dijet_boosted_mixed"}

        /*{1, "tt_0jet_mixed_idg0p5"},
        {2, "tt_boosted_mixed_idg0p5"},
        {3, "tt_dijet_mixed_idg0p5"},

        {4, "tt_0jet_mixed_!idg0p5"},
        {5, "tt_boosted_mixed_!idg0p5"},
        {6, "tt_dijet_mixed_!idg0p5"}*/
      };
      cats["tt_2016"] = {
        {31, "tt_0jet_other"},
        {32, "tt_boosted_other"},
        {33, "tt_dijet_lowboost_other"},
        {34, "tt_dijet_boosted_other"}
      };
    }
    else if (doDecays && do_mva && !id_cats) {
      cats_cp["tt_2016"] = {
        /*{1, "tt_higgs"},
        {2, "tt_zttEmbed"},
        {3, "tt_jetFakes"},*/
        //{4, "tt_misc"},
        //
        /*{1, "tt_higgs_rho"},
        {2, "tt_zttEmbed_rho"},
        {3, "tt_jetFakes_rho"},
        //
        {6, "tt_higgs_a1rho"},
        {7, "tt_zttEmbed_a1rho"},
        {8, "tt_jetFakes_a1rho"},*/

        // combined Higgs with MVA DM
        {1, "tt_higgs_mvarho"},
        {2, "tt_zttEmbed_mvarho"},
        {3, "tt_jetFakes_mvarho"},
        //
        {6, "tt_higgs_mvaa1rho"},
        {7, "tt_zttEmbed_mvaa1rho"},
        {8, "tt_jetFakes_mvaa1rho"},

        // split Higgs with MVA DM
        /*{1, "tt_ggh_mvarho"},
        {2, "tt_qqh_mvarho"},
        {3, "tt_zttEmbed_mvarho"},
        {4, "tt_jetFakes_mvarho"},
        //
        {6, "tt_ggh_mvaa1rho"},
        {7, "tt_qqh_mvaa1rho"},
        {8, "tt_zttEmbed_mvaa1rho"},
        {9, "tt_jetFakes_mvaa1rho"},*/

          // vienna NN
        /*{1, "tt_ggh_rho"},
        {2, "tt_qqh_rho"},
        {3, "tt_zttEmbed_rho"},
        {4, "tt_jetFakes_rho"},
        {5, "tt_misc_rho"},

        {6, "tt_ggh_a1rho"},
        {7, "tt_qqh_a1rho"},
        {8, "tt_zttEmbed_a1rho"},
        {9, "tt_jetFakes_a1rho"},
        {10, "tt_misc_a1rho"},*/

      };
      cats_cp["mt_2016"] = {
        {1,    "mt_higgs"}, 
        {2, "mt_zttEmbed"},
        {3, "mt_jetFakes"},
        {4, "mt_zll"},
        {5, "mt_tt"}
      };
      cats_cp["et_2016"] = {
        {1,    "et_higgs"}, 
        {2, "et_zttEmbed"},
        {3, "et_jetFakes"},
        {4, "et_zll"},
        {5, "et_tt"}
      };
      cats["tt_2016"] = {
        /*{31, "tt_higgs_other"},
        {32, "tt_zttEmbed_other"},
        {33, "tt_jetFakes_other"},*/
        //{34, "tt_misc_other"},
        //
        
        // with MVA DM
        {31, "tt_higgs_mvaother"},
        {32, "tt_zttEmbed_mvaother"},
        {33, "tt_jetFakes_mvaother"},

        // split Higgs with MVA DM
        /*{31, "tt_ggh_mvaother"},
        {32, "tt_qqh_mvaother"},
        {33, "tt_zttEmbed_mvaother"},
        {34, "tt_jetFakes_mvaother"},*/

          // vienna NN
        /*{31, "tt_ggh_other"},
        {32, "tt_qqh_other"},
        {33, "tt_zttEmbed_other"},
        {34, "tt_jetFakes_other"},
        {35, "tt_misc_other"},*/

      };
      cats["mt_2016"] = {
        {31,    "mt_higgs_other"}, 
        {32, "mt_zttEmbed_other"},
        {33, "mt_jetFakes_other"},
        {34, "mt_zll_other"},
        {35, "mt_tt_other"}
      };
      cats["et_2016"] = {
        {31,    "et_higgs_other"}, 
        {32, "et_zttEmbed_other"},
        {33, "et_jetFakes_other"},
        {34, "et_zll_other"},
        {35, "et_tt_other"}
      };
    }
    if (doDecays && do_mva && id_cats) {
      cats_cp["tt_2016"] = {
        {1, "tt_higgs"},
        {2, "tt_zttEmbed"},
        {3, "tt_jetFakes"}
      };
    }


    if(do_control_plots>0) {
      std::string extra="";
      if(do_control_plots==2) extra="lomsv_";
      if(do_control_plots==3) extra="himsv_";
      if(era=="2016"){
        cats_cp["et_2016"] = {};
        cats_cp["mt_2016"] = {};
        cats_cp["tt_2016"] = {};
        cats_cp["em_2016"] = {};
        cats["et_2016"] = {
          {100, "et_"+extra+"pt_1"},
          {101, "et_"+extra+"pt_2"},
          {102, "et_"+extra+"met"},
          {103, "et_"+extra+"pt_tt"},
          {104, "et_"+extra+"m_vis"},
          {105, "et_"+extra+"mjj"},
          {106, "et_"+extra+"sjdphi"},
          {107, "et_"+extra+"n_jets"}, 
          {108, "et_"+extra+"m_sv"}
         };
         cats["mt_2016"] = {
          {100, "mt_"+extra+"pt_1"},
          {101, "mt_"+extra+"pt_2"},
          {102, "mt_"+extra+"met"},
          {103, "mt_"+extra+"pt_tt"},
          {104, "mt_"+extra+"m_vis"},
          {105, "mt_"+extra+"mjj"},
          {106, "mt_"+extra+"sjdphi"},  
          {107, "mt_"+extra+"n_jets"},
          {108, "mt_"+extra+"m_sv"}
         };
         cats["tt_2016"] = {
          {100, "tt_"+extra+"n_jets"},
          /*{100, "tt_"+extra+"pt_1"},
          {101, "tt_"+extra+"pt_2"},
          {102, "tt_"+extra+"met"},
          {103, "tt_"+extra+"pt_tt"},
          {104, "tt_"+extra+"m_vis"},
          {105, "tt_"+extra+"mjj"},
          {106, "tt_"+extra+"sjdphi"},  
          {107, "tt_"+extra+"n_jets"},
          {108, "tt_"+extra+"m_sv"},*/
         };
         cats["em_2016"] = {
          {100, "em_"+extra+"pt_1"},
          {101, "em_"+extra+"pt_2"},
          {102, "em_"+extra+"met"},
          {103, "em_"+extra+"pt_tt"},
          {104, "em_"+extra+"m_vis"},
          {105, "em_"+extra+"mjj"},
          {106, "em_"+extra+"sjdphi"},
          {107, "em_"+extra+"n_jets"},
          {108, "em_"+extra+"m_sv"}
         }; 
       }
     if(era=="2017"){
        cats_cp["et_2017"] = {};
        cats_cp["mt_2017"] = {};
        cats_cp["tt_2017"] = {};
        cats_cp["em_2017"] = {};
        cats["et_2017"] = {
          {100, "et_"+extra+"pt_1"},
          {101, "et_"+extra+"pt_2"},
          {102, "et_"+extra+"met"},
          {103, "et_"+extra+"pt_tt"},
          {104, "et_"+extra+"m_vis"},
          {105, "et_"+extra+"mjj"},
          {106, "et_"+extra+"sjdphi"},
          {107, "et_"+extra+"n_jets"},
          {108, "et_"+extra+"m_sv"}
         };
         cats["mt_2017"] = {
          {100, "mt_"+extra+"pt_1"},
          {101, "mt_"+extra+"pt_2"},
          {102, "mt_"+extra+"met"},
          {103, "mt_"+extra+"pt_tt"},
          {104, "mt_"+extra+"m_vis"},
          {105, "mt_"+extra+"mjj"},
          {106, "mt_"+extra+"sjdphi"},
          {107, "mt_"+extra+"n_jets"},
          {108, "mt_"+extra+"m_sv"}
         };
         cats["tt_2017"] = {
          {100, "tt_"+extra+"pt_1"},
          {101, "tt_"+extra+"pt_2"},
          {102, "tt_"+extra+"met"},
          {103, "tt_"+extra+"pt_tt"},
          {104, "tt_"+extra+"m_vis"},
          {105, "tt_"+extra+"mjj"},
          {106, "tt_"+extra+"sjdphi"},
          {107, "tt_"+extra+"n_jets"},
          {108, "tt_"+extra+"m_sv"}
         };
         cats["em_2017"] = {
          {100, "em_"+extra+"pt_1"},
          {101, "em_"+extra+"pt_2"},
          {102, "em_"+extra+"met"},
          {103, "em_"+extra+"pt_tt"},
          {104, "em_"+extra+"m_vis"},
          {105, "em_"+extra+"mjj"},
          {106, "em_"+extra+"sjdphi"},
          {107, "em_"+extra+"n_jets"},
          {108, "em_"+extra+"m_sv"}
         };
       }

     }
    
    map<string, VString> sig_procs;
    sig_procs["ggH"] = {"ggH_sm_htt", "ggH_ps_htt", "ggH_mm_htt"};
    sig_procs["qqH"] = {"qqH_sm_htt", "qqH_ps_htt", "qqH_mm_htt"};   
 
    vector<string> masses = {"125"};    
    
    using ch::syst::bin_id;
    
    //! [part2]
    for(auto year: years) {
      for (auto chn : chns) {
          cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, cats[chn+"_"+year]);
          cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, cats_cp[chn+"_"+year]);

          cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, bkg_procs[chn], cats[chn+"_"+year], false);
          cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {chn+"_"+year}, bkg_procs[chn], cats_cp[chn+"_"+year], false);

          if(chn == "em" || chn == "et" || chn == "mt" || chn == "tt"){
            cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn+"_"+year}, sig_procs["qqH"], cats[chn+"_"+year], true); // SM VBF/VH are added as signal
            cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn+"_"+year}, sig_procs["qqH"], cats_cp[chn+"_"+year], true);

            cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn+"_"+year}, sig_procs["ggH"], cats[chn+"_"+year], true);
            cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn+"_"+year}, sig_procs["ggH"], cats_cp[chn+"_"+year], true);
          }
      }
    } 
    //! [part4]
    
    
    ch::AddSMRun2Systematics(cb, 0, ttbar_fit, false);
    
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
          if (year == "2017" && !do_control_plots) extra = "/2017/";
          if(chn == "ttbar") channel = "em"; 
          cb.cp().channel({chn+"_"+year}).backgrounds().ExtractShapes(
                                                             input_dir[chn] + extra + "htt_"+channel+".inputs-sm-13TeV"+postfix+".root",
                                                             "$BIN/$PROCESS",
                                                             "$BIN/$PROCESS_$SYSTEMATIC");
          if(chn == "em" || chn == "et" || chn == "mt" || chn == "tt"){
            cb.cp().channel({chn+"_"+year}).process(sig_procs["ggH"]).ExtractShapes(
                                                                    input_dir[chn] + extra + "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                                    "$BIN/$PROCESS$MASS",
                                                                    "$BIN/$PROCESS$MASS_$SYSTEMATIC");
            cb.cp().channel({chn+"_"+year}).process(sig_procs["qqH"]).ExtractShapes(
                                                                    input_dir[chn] + extra +  "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                                    "$BIN/$PROCESS$MASS",
                                                                    "$BIN/$PROCESS$MASS_$SYSTEMATIC");
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
    
    
    // And convert any shapes in the ttbar CRs to lnN:
    // Convert all shapes to lnN at this stage

    cb.cp().channel({"ttbar_2016","ttbar_2017"}).syst_type({"shape"}).ForEachSyst([](ch::Systematic *sys) {
        sys->set_type("lnN");
        if(sys->value_d() <0.001) {sys->set_value_d(0.001);};
        if(sys->value_u() <0.001) {sys->set_value_u(0.001);};
    });

    
    
    // can auto-merge the bins with bbb uncertainty > 90% - may be better to merge these bins by hand though!
    auto rebin = ch::AutoRebin()
    .SetBinThreshold(0.)
    .SetBinUncertFraction(0.9)
    .SetRebinMode(1)
    .SetPerformRebin(true)
    .SetVerbosity(1);
    if(auto_rebin) rebin.Rebin(cb, cb);
  
  
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


	
    // rename embedded energy-scale uncertainties so that they are not correlated with MC energy-scales
     cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_scale_e_13TeV","CMS_scale_embedded_e_13TeV"); 
     cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_scale_t_1prong_13TeV","CMS_scale_embedded_t_1prong_13TeV"); 
     cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_scale_t_1prong1pizero_13TeV","CMS_scale_embedded_t_1prong1pizero_13TeV"); 
     cb.cp().process({"EmbedZTT"}).RenameSystematic(cb,"CMS_scale_t_3prong_13TeV","CMS_scale_embedded_t_3prong_13TeV"); 


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
	
	
	if(do_control_plots==0) writer.WriteCards("cmb", cb);
	//Add all di-jet categories combined
	//
	
        writer.WriteCards("htt_2016", cb.cp().channel({"em_2016","et_2016","mt_2016","tt_2016","ttbar_2016"}));
        writer.WriteCards("htt_2017", cb.cp().channel({"em_2017","et_2017","mt_2017","tt_2017","ttbar_2017"})); 

	for (auto chn : cb.channel_set()) {
		 writer.WriteCards("htt_"+chn+"_dijet", cb.cp().channel({chn}).bin_id({3,4,5,6}));  

        // per-channel
	    writer.WriteCards(chn, cb.cp().channel({chn}));
        // And per-channel-category
        if (!do_mva) {
          writer.WriteCards("htt_"+chn+"_1_13TeV", cb.cp().channel({chn}).bin_id({1}));
          writer.WriteCards("htt_"+chn+"_2_13TeV", cb.cp().channel({chn}).bin_id({2}));
	      writer.WriteCards("htt_"+chn+"_3_13TeV", cb.cp().channel({chn}).bin_id({3}));
	      writer.WriteCards("htt_"+chn+"_4_13TeV", cb.cp().channel({chn}).bin_id({4}));
          writer.WriteCards("htt_"+chn+"_5_13TeV", cb.cp().channel({chn}).bin_id({5}));
          writer.WriteCards("htt_"+chn+"_6_13TeV", cb.cp().channel({chn}).bin_id({6}));
	}
        else {
          writer.WriteCards("htt_"+chn+"_1_13TeV", cb.cp().channel({chn}).bin_id({1}));
          writer.WriteCards("htt_"+chn+"_6_13TeV", cb.cp().channel({chn}).bin_id({6}));
          writer.WriteCards("htt_"+chn+"_7_13TeV", cb.cp().channel({chn}).bin_id({7}));
          writer.WriteCards("htt_"+chn+"_2_13TeV", cb.cp().channel({chn}).bin_id({2}));
          writer.WriteCards("htt_"+chn+"_3_13TeV", cb.cp().channel({chn}).bin_id({3}));
          writer.WriteCards("htt_"+chn+"_31_13TeV", cb.cp().channel({chn}).bin_id({31}));
          writer.WriteCards("htt_"+chn+"_32_13TeV", cb.cp().channel({chn}).bin_id({32}));
          writer.WriteCards("htt_"+chn+"_33_13TeV", cb.cp().channel({chn}).bin_id({33}));
          writer.WriteCards("htt_"+chn+"_34_13TeV", cb.cp().channel({chn}).bin_id({34}));
          writer.WriteCards("htt_"+chn+"_35_13TeV", cb.cp().channel({chn}).bin_id({35}));
          writer.WriteCards("htt_"+chn+"_36_13TeV", cb.cp().channel({chn}).bin_id({36}));
          writer.WriteCards("htt_"+chn+"_37_13TeV", cb.cp().channel({chn}).bin_id({37}));
          writer.WriteCards("htt_"+chn+"_38_13TeV", cb.cp().channel({chn}).bin_id({38}));
          writer.WriteCards("htt_"+chn+"_39_13TeV", cb.cp().channel({chn}).bin_id({39}));
          writer.WriteCards("htt_"+chn+"_41_13TeV", cb.cp().channel({chn}).bin_id({41}));
          writer.WriteCards("htt_"+chn+"_42_13TeV", cb.cp().channel({chn}).bin_id({42}));
          writer.WriteCards("htt_"+chn+"_43_13TeV", cb.cp().channel({chn}).bin_id({43}));
          writer.WriteCards("htt_"+chn+"_44_13TeV", cb.cp().channel({chn}).bin_id({44}));
          writer.WriteCards("htt_"+chn+"_45_13TeV", cb.cp().channel({chn}).bin_id({45}));
          writer.WriteCards("htt_"+chn+"_46_13TeV", cb.cp().channel({chn}).bin_id({46}));
          writer.WriteCards("htt_"+chn+"_47_13TeV", cb.cp().channel({chn}).bin_id({47}));
          writer.WriteCards("htt_"+chn+"_48_13TeV", cb.cp().channel({chn}).bin_id({48}));
          writer.WriteCards("htt_"+chn+"_49_13TeV", cb.cp().channel({chn}).bin_id({49}));

          writer.WriteCards("htt_"+chn+"_rho_13TeV", cb.cp().channel({chn}).bin_id({1,2,3,4,5}));
          writer.WriteCards("htt_"+chn+"_a1rho_13TeV", cb.cp().channel({chn}).bin_id({6,7,8,9,10}));
          writer.WriteCards("htt_"+chn+"_others_13TeV", cb.cp().channel({chn}).bin_id({31,32,33,34}));

        }
        
    }
    
    
    cb.PrintAll();
    cout << " done\n";
    
    
}
