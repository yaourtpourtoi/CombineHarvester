{

    TString path="output/tauIDSF_output/MC/2018/htt_mt_1_13TeV/125/multidimfit.r.root";
    
    TFile *f = new TFile(path, "read");
    
    RooArgList result_list=fit_mdf->floatParsFinal();
    
    TH1F *hist = new TH1F("hist","hist", result_list.getSize()-2, 0 ,result_list.getSize()-2);
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();
    for (int i=0; i<result_list.getSize()-2; i++){
        RooRealVar *result = (RooRealVar*) result_list.at(i);
        //cout<<result->GetName()<<"--------" <<result->getValV()<<"--------"<<result->getError()<<endl;
        //hist->SetBinContent(i+1, result->getValV());
        hist->SetBinContent(i+1, 10);
        hist->SetBinError(i+1, result->getError());
        hist->GetXaxis()->SetBinLabel(i+1, result->GetName());
        
        gr->SetPoint(i+1, result->getValV(), i+0.5);
        gr->SetPointError(i+1, result->getError(), result->getError(), 0, 0);
        //gr->GetXaxis()->SetBinLabel(i+1, result->GetName());
    }
        hist->SetBarWidth(0);
        hist->SetMaximum(3);
        hist->SetMinimum(-3);
        hist->SetStats(0);
        hist->SetFillStyle(3359);
        hist->GetYaxis()->SetTitle("(#theta-#theta_{0})/#sigma_{pre-fit}");
        
    
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(kBlue);
    gr->SetLineWidth(2);
    gr->SetLineColor(kBlue);
 TCanvas *c1 = new TCanvas("c1","c1",800,800);
    gPad->SetLeftMargin(0.3);
 hist->Draw("E2 hbar");
    gr->Draw("EPsame");
    c1->SetGridx();
    //list_result = fit_mdf->floatParsFinal();
    
    //cout<<list_result.getSize();



}
