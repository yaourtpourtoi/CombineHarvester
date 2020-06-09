{

//-----------Only this part should be modified
int year =2016;//2016 or 2017 or  2018
TString sample = "embed"; //MC or embed
TString channel = "et"; // "ttAndmt" or "et"
//--------------------------------

gROOT->SetBatch(1);

    const Int_t n = 18; 

    for (int dm=0; dm<n; dm++){
        
        stringstream fileNameStream;
        fileNameStream <<"output/"<<channel<<"_datacards_output/"<<sample<<"/"<<to_string(year)<<"/htt_mt_"<<dm+1<<"_13TeV/125/multidimfit.r.root"<<endl;

        TString path;
        fileNameStream >> path;

    
    TFile *f = new TFile(path, "read");
    
    RooArgList result_list=fit_mdf->floatParsFinal();
    
    TH1F *hist = new TH1F("hist","", result_list.getSize()-2, 0 ,result_list.getSize()-2);
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();
    for (int i=0; i<result_list.getSize()-2; i++){
        RooRealVar *result = (RooRealVar*) result_list.at(i);
        hist->SetBinContent(i+1, 10);
        hist->SetBinError(i+1, result->getError());
        hist->GetXaxis()->SetBinLabel(i+1, result->GetName());
        
        gr->SetPoint(i+1, result->getValV(), i+0.5);
        gr->SetPointError(i+1, result->getError(), result->getError(), 0, 0);
    }
        hist->SetBarWidth(0);
        hist->SetMaximum(3);
        hist->SetMinimum(-3);
        hist->SetStats(0);
        //hist->SetFillStyle(3359);
        hist->GetYaxis()->SetTitle("(#theta-#theta_{0})/#sigma_{pre-fit}");
        
   hist->GetXaxis()->SetLabelSize(0.03); 
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(kBlue);
    gr->SetLineWidth(2);
    gr->SetLineColor(kBlue);
    gr->SetTitle("");
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
//    c1->SetBatch(kTrue);
    gPad->SetLeftMargin(0.3);
    hist->Draw("E2 hbar");
    gr->Draw("EPsame");
    c1->SetGridx();

stringstream fileOutStream;
fileOutStream<<"output_pulls/"<<channel<<"/"<<sample<<"/"<<to_string(year)<<"/pull_"<<dm+1<<".pdf"<<endl;
TString out_dir;
fileOutStream>>out_dir;

    c1->SaveAs(out_dir);
    

}

}
