{
   TFile *file1 = new TFile("tauIDSFGraphs.root","RECREATE");
gROOT->SetBatch();
   auto c1 = new TCanvas("c1","c1",200,10,700,500);
   
   
   
   
   const Int_t n = 5;


   Double_t x[n] = {1, 2, 3, 4, 5};
   Double_t ex[n] = {0, 0, 0, 0, 0};


//embedlowpt
   Double_t yEmbedLow[n]  = {1.0342, 9.5716E-01, 1.0011, 9.3268E-01, 8.8146E-01};
   Double_t eyEmbedLow[n]  = {1.95E-02, 1.31E-02, 2.52E-02, 1.89E-02, 3.57E-02};
   
//embedhighpt
   Double_t yEmbedHigh[n]  = {1.0584, 9.2757E-01, 9.941E-01, 9.577E-01, 8.9251E-01};
   Double_t eyEmbedHigh[n]  = {4.36E-02, 2.63E-02, 4.83E-02, 3.27E-02, 4.98E-02};

//MClowpt
   Double_t yMCLow[n]  = {9.7561E-01, 9.4249E-01, 8.9611E-01, 9.01E-01, 8.4193E-01};
  Double_t eyMCLow[n]  = {6.41E-03, 9.91E-03, 1.85E-02, 1.41E-02, 2.62E-02};

//MChighpt
   Double_t yMCHigh[n]  = {1.0136, 9.7727E-01, 9.5469E-01, 9.4657E-01, 9.1228E-01};
   Double_t eyMCHigh[n]  = {3.36E-02, 2.17E-02, 3.19E-02, 2.31E-02, 3.99E-02};
   
   auto gr_embed_20pt40 = new TGraphErrors(n,x,yEmbedLow,ex,eyEmbedLow);
   gr_embed_20pt40->GetHistogram()->SetMaximum(1.2);
   gr_embed_20pt40->GetHistogram()->SetMinimum(0.8);
   gr_embed_20pt40->SetTitle("Embedding (20<pt<40)");
   gr_embed_20pt40->GetYaxis()->SetTitle("Scale Factor");
   gr_embed_20pt40->GetXaxis()->SetBinLabel(gr_embed_20pt40->GetXaxis()->FindBin(1),"DM0");
   gr_embed_20pt40->GetXaxis()->SetBinLabel(gr_embed_20pt40->GetXaxis()->FindBin(2),"DM1");
   gr_embed_20pt40->GetXaxis()->SetBinLabel(gr_embed_20pt40->GetXaxis()->FindBin(3),"DM2");
   gr_embed_20pt40->GetXaxis()->SetBinLabel(gr_embed_20pt40->GetXaxis()->FindBin(4),"DM10");
   gr_embed_20pt40->GetXaxis()->SetBinLabel(gr_embed_20pt40->GetXaxis()->FindBin(5),"DM11");
   gr_embed_20pt40->GetXaxis()->LabelsOption("h");
   gr_embed_20pt40->GetXaxis()->SetLabelSize(0.07);
   gr_embed_20pt40->GetYaxis()->SetTitleSize(0.05);
   gr_embed_20pt40->GetYaxis()->SetTitleOffset(0.8);
   gr_embed_20pt40->SetMarkerColor(4);
   gr_embed_20pt40->SetMarkerStyle(21);
   gr_embed_20pt40->Draw("AP");

gr_embed_20pt40->Write("gr_embed_20pt40");




   auto c2 = new TCanvas("c2","c2",200,10,700,500);
   auto gr_embed_40pt = new TGraphErrors(n,x,yEmbedHigh,ex,eyEmbedHigh);
   gr_embed_40pt->GetHistogram()->SetMaximum(1.2);
   gr_embed_40pt->GetHistogram()->SetMinimum(0.8);
   gr_embed_40pt->SetTitle("Embedding (40<pt)");
   gr_embed_40pt->GetYaxis()->SetTitle("Scale Factor");
   gr_embed_40pt->GetXaxis()->SetBinLabel(gr_embed_40pt->GetXaxis()->FindBin(1),"DM0");
   gr_embed_40pt->GetXaxis()->SetBinLabel(gr_embed_40pt->GetXaxis()->FindBin(2),"DM1");
   gr_embed_40pt->GetXaxis()->SetBinLabel(gr_embed_40pt->GetXaxis()->FindBin(3),"DM2");
   gr_embed_40pt->GetXaxis()->SetBinLabel(gr_embed_40pt->GetXaxis()->FindBin(4),"DM10");
   gr_embed_40pt->GetXaxis()->SetBinLabel(gr_embed_40pt->GetXaxis()->FindBin(5),"DM11");
   gr_embed_40pt->GetXaxis()->LabelsOption("h");
   gr_embed_40pt->GetXaxis()->SetLabelSize(0.07);
   gr_embed_40pt->GetYaxis()->SetTitleSize(0.05);
   gr_embed_40pt->GetYaxis()->SetTitleOffset(0.8);
   gr_embed_40pt->SetMarkerColor(4);
   gr_embed_40pt->SetMarkerStyle(21);
   gr_embed_40pt->Draw("AP");

gr_embed_40pt->Write("gr_embed_40pt");

   auto c3 = new TCanvas("c3","c3",200,10,700,500);
   auto gr_MC_20pt40 = new TGraphErrors(n,x,yMCLow,ex,eyMCLow);
   gr_MC_20pt40->GetHistogram()->SetMaximum(1.2);
   gr_MC_20pt40->GetHistogram()->SetMinimum(0.8);
   gr_MC_20pt40->SetTitle("MC (20<pt<40)");
   gr_MC_20pt40->GetYaxis()->SetTitle("Scale Factor");
   gr_MC_20pt40->GetXaxis()->SetBinLabel(gr_MC_20pt40->GetXaxis()->FindBin(1),"DM0");
   gr_MC_20pt40->GetXaxis()->SetBinLabel(gr_MC_20pt40->GetXaxis()->FindBin(2),"DM1");
   gr_MC_20pt40->GetXaxis()->SetBinLabel(gr_MC_20pt40->GetXaxis()->FindBin(3),"DM2");
   gr_MC_20pt40->GetXaxis()->SetBinLabel(gr_MC_20pt40->GetXaxis()->FindBin(4),"DM10");
   gr_MC_20pt40->GetXaxis()->SetBinLabel(gr_MC_20pt40->GetXaxis()->FindBin(5),"DM11");
   gr_MC_20pt40->GetXaxis()->LabelsOption("h");
   gr_MC_20pt40->GetXaxis()->SetLabelSize(0.07);
   gr_MC_20pt40->GetYaxis()->SetTitleSize(0.05);
   gr_MC_20pt40->GetYaxis()->SetTitleOffset(0.8);
   gr_MC_20pt40->SetMarkerColor(4);
   gr_MC_20pt40->SetMarkerStyle(21);
   gr_MC_20pt40->Draw("AP");
 gr_MC_20pt40->Write("gr_MC_20pt40"); 
  
   auto c4 = new TCanvas("c4","c4",200,10,700,500);
   auto gr_MC_40pt = new TGraphErrors(n,x,yMCHigh,ex,eyMCHigh);
   gr_MC_40pt->GetHistogram()->SetMaximum(1.2);
   gr_MC_40pt->GetHistogram()->SetMinimum(0.8);
   gr_MC_40pt->SetTitle("MC (40<pt)");
   gr_MC_40pt->GetYaxis()->SetTitle("Scale Factor");
   gr_MC_40pt->GetXaxis()->SetBinLabel(gr_MC_40pt->GetXaxis()->FindBin(1),"DM0");
   gr_MC_40pt->GetXaxis()->SetBinLabel(gr_MC_40pt->GetXaxis()->FindBin(2),"DM1");
   gr_MC_40pt->GetXaxis()->SetBinLabel(gr_MC_40pt->GetXaxis()->FindBin(3),"DM2");
   gr_MC_40pt->GetXaxis()->SetBinLabel(gr_MC_40pt->GetXaxis()->FindBin(4),"DM10");
   gr_MC_40pt->GetXaxis()->SetBinLabel(gr_MC_40pt->GetXaxis()->FindBin(5),"DM11");
   gr_MC_40pt->GetXaxis()->LabelsOption("h");
   gr_MC_40pt->GetXaxis()->SetLabelSize(0.07);
   gr_MC_40pt->GetYaxis()->SetTitleSize(0.05);
   gr_MC_40pt->GetYaxis()->SetTitleOffset(0.8);
   gr_MC_40pt->SetMarkerColor(4);
   gr_MC_40pt->SetMarkerStyle(21);
   gr_MC_40pt->Draw("AP");

gr_MC_40pt->Write("gr_MC_40pt");




file1->Close();



//old numbers with smaller m_vis range for mt channel when producing datacards:
//
//embedlowpt
//   Double_t y[n]  = {1.0257, 0.944, 1.0035, 9.2362E-01, 8.7715E-01};
//   Double_t ey[n]  = {2.8E-02, 1.9E-02, 3.13E-02, 2.66E-02, 4.2E-02};
   
//embedhighpt
 //  Double_t y[n]  = {1.0571, 9.1499E-01, 9.9081E-01, 9.5522E-01, 8.8818E-01};
 //  Double_t ey[n]  = {4.33E-02, 2.79E-02, 4.64E-02, 3.55E-02, 5.24E-02};

//MClowpt
 //  Double_t y[n]  = {9.5595E-01,9.3059E-01,8.988E-01,8.9201E-01, 8.3738E-01};
  // Double_t ey[n]  = {1.63E-02, 9.91E-03, 1.88E-02, 1.48E-02, 2.94E-02};

//MChighpt
 //  Double_t y[n]  = {1.0209, 9.7049E-01, 9.5921E-01, 9.4256E-01, 9.1849E-01};
  // Double_t ey[n]  = {3.15E-02, 2.02E-02, 3.05E-02, 2.24E-02, 3.99E-02};

}

