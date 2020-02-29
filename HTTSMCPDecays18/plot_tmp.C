{
   auto c1 = new TCanvas("c1","c1",200,10,700,500);
   //c1->SetFillColor(42);
   //c1->SetGrid();
   //c1->GetFrame()->SetFillColor(21);
   //c1->GetFrame()->SetBorderSize(12);
   const Int_t n = 5;


   Double_t x[n] = {1, 2, 3, 4, 5};
   Double_t ex[n] = {0, 0, 0, 0, 0};


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
   Double_t y[n]  = {1.0209, 9.7049E-01, 9.5921E-01, 9.4256E-01, 9.1849E-01};
   Double_t ey[n]  = {3.15E-02, 2.02E-02, 3.05E-02, 2.24E-02, 3.99E-02};
   
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->GetHistogram()->SetMaximum(1.2);
   gr->GetHistogram()->SetMinimum(0.8);
   //gr->SetTitle("Embedding (20<pt<40)");
  // gr->SetTitle("Embedding (40<pt)");
  // gr->SetTitle("MC (20<pt<40)");
   gr->SetTitle("MC (40<pt)");
   gr->GetYaxis()->SetTitle("Scale Factor");
   gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(1),"DM0");
   gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(2),"DM1");
   gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(3),"DM2");
   gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(4),"DM10");
   gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(5),"DM11");
   gr->GetXaxis()->LabelsOption("h");
   gr->GetXaxis()->SetLabelSize(0.07);
   gr->GetYaxis()->SetTitleSize(0.05);
   gr->GetYaxis()->SetTitleOffset(0.8);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
   //c1->Print("TauIDSF_Embed_20pt40.pdf");
  /// c1->Print("TauIDSF_Embed_40pt.pdf");
  // c1->Print("TauIDSF_MC_20pt40.pdf");
   c1->Print("TauIDSF_MC_40pt.pdf");
}

