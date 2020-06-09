{

//-----------Only this part should be modified
TString sample = "MC"; //MC or embed
TString channel = "et"; // "ttAndmt" or "et"
//--------------------------------

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

int year_max=2018;
int year_min=2016;


Double_t ex[n] = {0, 0, 0, 0, 0};

Double_t SF_lowpt[n] = {0, 0,0, 0, 0};
Double_t eSF_lowpt[n] = {0, 0, 0, 0, 0};

Double_t SF_highpt[n] = {0, 0,0, 0, 0};
Double_t eSF_highpt[n] = {0, 0, 0, 0, 0};

TGraphErrors *gr_2018_lowpt = new TGraphErrors();
TGraphErrors *gr_2017_lowpt = new TGraphErrors();
TGraphErrors *gr_2016_lowpt = new TGraphErrors();

TGraphErrors *gr_2018_highpt = new TGraphErrors();
TGraphErrors *gr_2017_highpt = new TGraphErrors();
TGraphErrors *gr_2016_highpt = new TGraphErrors();

for (int year = year_max; year>=year_min; year--){
    for (int i=0; i<n; i++){
        
        stringstream fileNameStream_lowpt;
        fileNameStream_lowpt <<"output/"<<channel<<"_datacards_output/"<<sample<<"/"<<to_string(year)<<"/htt_mt_"<<i+1<<"_13TeV/125/multidimfit.r.root"<<endl;
        TString fileName_lowpt;
        fileNameStream_lowpt >> fileName_lowpt;
        
        stringstream fileNameStream_highpt;
        fileNameStream_highpt <<"output/"<<channel<<"_datacards_output/"<<sample<<"/"<<to_string(year)<<"/htt_mt_"<<5+i+1<<"_13TeV/125/multidimfit.r.root"<<endl;
        TString fileName_highpt;
        fileNameStream_highpt >> fileName_highpt;
        
        
        TFile *f_lowpt= new TFile(fileName_lowpt, "read");
        RooRealVar *r_result_lowpt = (RooRealVar*) fit_mdf->floatParsFinal().find("r");
        SF_lowpt[i] = r_result_lowpt->getValV();
        eSF_lowpt[i] = r_result_lowpt->getError();
        
        TFile *f_highpt= new TFile(fileName_highpt, "read");
        RooRealVar *r_result_highpt = (RooRealVar*) fit_mdf->floatParsFinal().find("r");
        SF_highpt[i] = r_result_highpt->getValV();
        eSF_highpt[i] = r_result_highpt->getError();
    }
    
    if (year==2018){
        gr_2018_lowpt = new TGraphErrors(n,x2,SF_lowpt,ex,eSF_lowpt);
        gr_2018_highpt = new TGraphErrors(n,x2,SF_highpt,ex,eSF_highpt);
    }
    if (year==2017){
        gr_2017_lowpt = new TGraphErrors(n,x1,SF_lowpt,ex,eSF_lowpt);
        gr_2017_highpt = new TGraphErrors(n,x1,SF_highpt,ex,eSF_highpt);
    }
    if (year==2016){
        gr_2016_lowpt = new TGraphErrors(n,x3,SF_lowpt,ex,eSF_lowpt);
        gr_2016_highpt = new TGraphErrors(n,x3,SF_highpt,ex,eSF_highpt);
    }
    
}

max_gr =1.2;
min_gr = 0.7;

   TCanvas *c1 = new TCanvas();
   gr_2018_lowpt->GetHistogram()->SetMaximum(max_gr);
   gr_2018_lowpt->GetHistogram()->SetMinimum(min_gr);
   gr_2018_lowpt->GetYaxis()->SetTitle("Scale Factor");
   gr_2018_lowpt->GetXaxis()->SetBinLabel(gr_2018_lowpt->GetXaxis()->FindBin(1),"DM0");
   gr_2018_lowpt->GetXaxis()->SetBinLabel(gr_2018_lowpt->GetXaxis()->FindBin(2),"DM1");
   gr_2018_lowpt->GetXaxis()->SetBinLabel(gr_2018_lowpt->GetXaxis()->FindBin(3),"DM2");
   gr_2018_lowpt->GetXaxis()->SetBinLabel(gr_2018_lowpt->GetXaxis()->FindBin(4),"DM10");
   gr_2018_lowpt->GetXaxis()->SetBinLabel(gr_2018_lowpt->GetXaxis()->FindBin(5),"DM11");
   gr_2018_lowpt->GetXaxis()->LabelsOption("h");
   gr_2018_lowpt->GetXaxis()->SetLabelSize(0.07);
   gr_2018_lowpt->GetYaxis()->SetTitleSize(0.05);
   gr_2018_lowpt->GetYaxis()->SetTitleOffset(0.8);
   
   gr_2018_lowpt->SetMarkerColor(kBlue);
   gr_2018_lowpt->SetMarkerStyle(20);
   gr_2018_lowpt->SetLineColor(kBlue);

   gr_2017_lowpt->SetMarkerColor(kRed);
   gr_2017_lowpt->SetMarkerStyle(20);
   gr_2017_lowpt->SetLineColor(kRed);
 
   gr_2016_lowpt->SetMarkerColor(kBlack);
   gr_2016_lowpt->SetMarkerStyle(20);
   gr_2016_lowpt->SetLineColor(kBlack);

auto mg_lowpt = new TMultiGraph();

mg_lowpt->Add(gr_2017_lowpt);
mg_lowpt->Add(gr_2016_lowpt);


leg1= new TLegend(0.55,0.65,0.7,0.85);
leg1->AddEntry(gr_2018_lowpt,"2018","lep");
leg1->AddEntry(gr_2017_lowpt,"2017","lep");
leg1->AddEntry(gr_2016_lowpt,"2016","lep");
leg1->SetFillColor(kWhite);



mg_lowpt->SetTitle("");
gr_2018_lowpt->SetTitle("");

gr_2018_lowpt->Draw("AP");
mg_lowpt->Draw("P");

leg1->Draw();

stringstream outputstream_lowpt;
outputstream_lowpt<<"output_TauIDresult_Generic/"<<channel<<"_"<<sample<<"_20pt40_MVADM.pdf"<<endl;

TString output_lowpt;
outputstream_lowpt>>output_lowpt;

c1->Print(output_lowpt);



   TCanvas *c2 = new TCanvas();
   gr_2018_highpt->GetHistogram()->SetMaximum(max_gr);
   gr_2018_highpt->GetHistogram()->SetMinimum(min_gr);
   gr_2018_highpt->GetYaxis()->SetTitle("Scale Factor");
   gr_2018_highpt->GetXaxis()->SetBinLabel(gr_2018_highpt->GetXaxis()->FindBin(1),"DM0");
   gr_2018_highpt->GetXaxis()->SetBinLabel(gr_2018_highpt->GetXaxis()->FindBin(2),"DM1");
   gr_2018_highpt->GetXaxis()->SetBinLabel(gr_2018_highpt->GetXaxis()->FindBin(3),"DM2");
   gr_2018_highpt->GetXaxis()->SetBinLabel(gr_2018_highpt->GetXaxis()->FindBin(4),"DM10");
   gr_2018_highpt->GetXaxis()->SetBinLabel(gr_2018_highpt->GetXaxis()->FindBin(5),"DM11");
   gr_2018_highpt->GetXaxis()->LabelsOption("h");
   gr_2018_highpt->GetXaxis()->SetLabelSize(0.07);
   gr_2018_highpt->GetYaxis()->SetTitleSize(0.05);
   gr_2018_highpt->GetYaxis()->SetTitleOffset(0.8);
   
   gr_2018_highpt->SetMarkerColor(kBlue);
   gr_2018_highpt->SetMarkerStyle(20);
   gr_2018_highpt->SetLineColor(kBlue);

   gr_2017_highpt->SetMarkerColor(kRed);
   gr_2017_highpt->SetMarkerStyle(20);
   gr_2017_highpt->SetLineColor(kRed);
 
   gr_2016_highpt->SetMarkerColor(kBlack);
   gr_2016_highpt->SetMarkerStyle(20);
   gr_2016_highpt->SetLineColor(kBlack);

auto mg_highpt = new TMultiGraph();

mg_highpt->Add(gr_2017_highpt);
mg_highpt->Add(gr_2016_highpt);


leg2= new TLegend(0.55,0.65,0.7,0.85);
leg2->AddEntry(gr_2018_highpt,"2018","lep");
leg2->AddEntry(gr_2017_highpt,"2017","lep");
leg2->AddEntry(gr_2016_highpt,"2016","lep");
leg2->SetFillColor(kWhite);



mg_highpt->SetTitle("");
gr_2018_highpt->SetTitle("");

gr_2018_highpt->Draw("AP");
mg_highpt->Draw("P");

leg2->Draw();



stringstream outputstream_highpt;
outputstream_highpt<<"output_TauIDresult_Generic/"<<channel<<"_"<<sample<<"_40pt_MVADM.pdf"<<endl;

TString output_highpt;
outputstream_highpt>>output_highpt;

c2->Print(output_highpt);



}





