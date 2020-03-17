{

//path = "output/tauIDSF_output/MC/2017/htt_mt_9_13TeV/125/multidimfit.r.root";
//
//
//path1= "output/tauIDSF_output/MC/2017/htt_mt_1_13TeV/125/multidimfit.r.root";
//
//
//string i="1";
//fileNameStream <<"output/tauIDSF_output/MC/2017/htt_mt_"<<i<<"_13TeV/125/multidimfit.r.root"<<endl;
//TString fileName;
//fileNameStream >> fileName;
//cout<<fileName<<endl;

    //stringstream fileNameStream;
    //fileNameStream <<"output/tauIDSF_output/MC/2017/htt_mt_"<<1<<"_13TeV/125/multidimfit.r.root"<<endl;
    //TString fileName;
    //fileNameStream >> fileName;
    //TFile f(fileName, "read");
    //result = f.Get("fit_mdf");


const Int_t n = 5; 

Double_t x1[n] = {1, 2, 3, 4, 5};
Double_t x2[n] = {1, 2, 3, 4, 5};
Double_t x3[n] = {1, 2, 3, 4, 5};

offset=0.1;

for (int k=0;k<n;k++)
{
x2[k]+=offset;
x3[k]-=offset;
}



Double_t ex[n] = {0, 0, 0, 0, 0};

Double_t SF[n] = {0, 0,0, 0, 0};
Double_t eSF[n] = {0, 0, 0, 0, 0};


int  year =2017;



TGraphErrors *gr_2018 = new TGraphErrors();
TGraphErrors *gr_2017 = new TGraphErrors();
TGraphErrors *gr_2016 = new TGraphErrors();


for (int year = 2018; year>2015; year--){
    for (int i=0; i<n; i++){
        cout<<"i= "<<i<<endl;
        
        stringstream fileNameStream;
        //pt>40:
        fileNameStream <<"output/tauIDSF_output/embed/"<<to_string(year)<<"/htt_mt_"<<5+i+1<<"_13TeV/125/multidimfit.r.root"<<endl;
        
        //20<pt<40:
        //fileNameStream <<"output/tauIDSF_output/MC/"<<to_string(year)<<"/htt_mt_"<<i+1<<"_13TeV/125/multidimfit.r.root"<<endl;
        
        TString fileName;
        fileNameStream >> fileName;
        
        cout<<fileName<<endl;
        TFile *f= new TFile(fileName, "read");
        //result = f->Get("fit_mdf");
        
        RooRealVar *r_result = (RooRealVar*) fit_mdf->floatParsFinal().find("r");
        
        SF[i] = r_result->getValV();
        eSF[i] = r_result->getError();
    }
    if (year==2018)
    gr_2018 = new TGraphErrors(n,x2,SF,ex,eSF);
    if (year==2017)
    gr_2017 = new TGraphErrors(n,x1,SF,ex,eSF);
    if (year==2016)
    gr_2016 = new TGraphErrors(n,x3,SF,ex,eSF);

}


gr_2018->Print();
cout<<"-------"<<endl;
gr_2017->Print();
cout<<"-------"<<endl;
gr_2016->Print();
cout<<"-------"<<endl;


//for (int j=0;j<n;j++)
//{cout<<SF[j]<<endl;}

max_gr =1.2;
min_gr = 0.7;

   TCanvas *c1 = new TCanvas();
   gr_2018->GetHistogram()->SetMaximum(max_gr);
   gr_2018->GetHistogram()->SetMinimum(min_gr);
   gr_2018->GetYaxis()->SetTitle("Scale Factor");
   gr_2018->GetXaxis()->SetBinLabel(gr_2018->GetXaxis()->FindBin(1),"DM0");
   gr_2018->GetXaxis()->SetBinLabel(gr_2018->GetXaxis()->FindBin(2),"DM1");
   gr_2018->GetXaxis()->SetBinLabel(gr_2018->GetXaxis()->FindBin(3),"DM2");
   gr_2018->GetXaxis()->SetBinLabel(gr_2018->GetXaxis()->FindBin(4),"DM10");
   gr_2018->GetXaxis()->SetBinLabel(gr_2018->GetXaxis()->FindBin(5),"DM11");
   gr_2018->GetXaxis()->LabelsOption("h");
   gr_2018->GetXaxis()->SetLabelSize(0.07);
   gr_2018->GetYaxis()->SetTitleSize(0.05);
   gr_2018->GetYaxis()->SetTitleOffset(0.8);
   
   
   
   gr_2018->SetMarkerColor(kBlue);
   gr_2018->SetMarkerStyle(20);
   gr_2018->SetLineColor(kBlue);

   gr_2017->SetMarkerColor(kRed);
   gr_2017->SetMarkerStyle(20);
   gr_2017->SetLineColor(kRed);
   
   gr_2016->SetMarkerColor(kBlack);
   gr_2016->SetMarkerStyle(20);
   gr_2016->SetLineColor(kBlack);


auto mg = new TMultiGraph();

//mg->Add(gr_2018);
mg->Add(gr_2017);
mg->Add(gr_2016);


leg= new TLegend(0.55,0.65,0.7,0.85);
leg->AddEntry(gr_2016,"2016","lep");
leg->AddEntry(gr_2017,"2017","lep");
leg->AddEntry(gr_2018,"2018","lep");
leg->SetFillColor(kWhite);



mg->SetTitle("");
gr_2018->SetTitle("MC - 40<pt - MVA DM");

gr_2018->Draw("AP");
mg->Draw("P");

leg->Draw();



c1->Print("output_result/VVLoose_MC_pt40_MVADM.pdf");





}





