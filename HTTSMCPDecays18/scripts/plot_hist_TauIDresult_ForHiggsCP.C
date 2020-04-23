{
//Specified for H->\tau\tau CP analysis: all HPS DM=0 events are excluded from MVA DM 1 and MVA DM 2


//-----------Only this part should be modified
int year =2016;//2016 or 2017 or  2018
TString sample = "embed"; //MC or embed
TString channel = "et"; // "ttAndmt" or "et"
//--------------------------------





Int_t n =22;    
double SF[22]={};
double eSF[22]={};

    for (int i=0; i<n; i++){
        
        stringstream fileNameStream;
        
        fileNameStream <<"output/"<<channel<<"_datacards_output/"<<sample<<"/"<<to_string(year)<<"/htt_mt_"<<i+1<<"_13TeV/125/multidimfit.r.root"<<endl;
        
        TString fileName;
        fileNameStream >> fileName;
        
        cout<<fileName<<endl;
        TFile *f= new TFile(fileName, "read");
        
        RooRealVar *r_result = (RooRealVar*) fit_mdf->floatParsFinal().find("r");
        
        SF[i] = r_result->getValV();
        eSF[i] = r_result->getError();
    }



TH1F *h_MVA_lowpt = new TH1F("h_MVA_lowpt","h_MVA_lowpt", 12, 0, 12);

    h_MVA_lowpt->SetBinContent(1,  SF[0]);
    h_MVA_lowpt->SetBinContent(2,  SF[18]);//MVADM==1 and HPSDM!=0
    h_MVA_lowpt->SetBinContent(3,  SF[19]);//MVADM==2 and HPSDM!=0
    h_MVA_lowpt->SetBinContent(11, SF[3]);
    h_MVA_lowpt->SetBinContent(12, SF[4]);    
    
    h_MVA_lowpt->SetBinError(1,  eSF[0]);
    h_MVA_lowpt->SetBinError(2,  eSF[18]);
    h_MVA_lowpt->SetBinError(3,  eSF[19]);
    h_MVA_lowpt->SetBinError(11, eSF[3]);
    h_MVA_lowpt->SetBinError(12, eSF[4]);    


TH1F *h_MVA_highpt = new TH1F("h_MVA_highpt","h_MVA_highpt", 12, 0, 12);

    h_MVA_highpt->SetBinContent(1,  SF[5]);
    h_MVA_highpt->SetBinContent(2,  SF[20]);
    h_MVA_highpt->SetBinContent(3,  SF[21]);
    h_MVA_highpt->SetBinContent(11, SF[8]);
    h_MVA_highpt->SetBinContent(12, SF[9]);    
    
    h_MVA_highpt->SetBinError(1,  eSF[5]);
    h_MVA_highpt->SetBinError(2,  eSF[20]);
    h_MVA_highpt->SetBinError(3,  eSF[21]);
    h_MVA_highpt->SetBinError(11, eSF[8]);
    h_MVA_highpt->SetBinError(12, eSF[9]);    


TH1F *h_HPS_lowpt = new TH1F("h_HPS_lowpt","h_HPS_lowpt", 12, 0, 12);

    h_HPS_lowpt->SetBinContent(1,  SF[10]);
    h_HPS_lowpt->SetBinContent(2,  SF[11]);
    h_HPS_lowpt->SetBinContent(11, SF[12]);
    h_HPS_lowpt->SetBinContent(12, SF[13]);
    

    h_HPS_lowpt->SetBinError(1,  eSF[10]);
    h_HPS_lowpt->SetBinError(2,  eSF[11]);
    h_HPS_lowpt->SetBinError(11, eSF[12]);
    h_HPS_lowpt->SetBinError(12, eSF[13]);


TH1F *h_HPS_highpt = new TH1F("h_HPS_highpt","h_HPS_highpt", 12, 0, 12);

    h_HPS_highpt->SetBinContent(1,  SF[14]);
    h_HPS_highpt->SetBinContent(2,  SF[15]);
    h_HPS_highpt->SetBinContent(11, SF[16]);
    h_HPS_highpt->SetBinContent(12, SF[17]);
    

    h_HPS_highpt->SetBinError(1,  eSF[14]);
    h_HPS_highpt->SetBinError(2,  eSF[15]);
    h_HPS_highpt->SetBinError(11, eSF[16]);
    h_HPS_highpt->SetBinError(12, eSF[17]);

stringstream rootFileStream;
rootFileStream <<"output_TauIDresult_ForHiggsCP/result_TauIDSF_"<<channel<<"_"<<sample<<"_"<<to_string(year)<<".root"<<endl;
TString rootFile;
rootFileStream >> rootFile;
cout<<rootFile<<endl;

TFile *file = new TFile(rootFile,"RECREATE");

h_MVA_lowpt->Write();
h_MVA_highpt->Write();
h_HPS_lowpt->Write();
h_HPS_highpt->Write();

file->Close();

}
