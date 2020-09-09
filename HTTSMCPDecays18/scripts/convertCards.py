import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

cp_bins = {

        "tt_201$_zttEmbed": 1,
        "tt_201$_jetFakes": 1,
        "tt_201$_higgs_Rho_Rho": 10,
        "tt_201$_higgs_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_higgs_A1_Rho": 4,
        "tt_201$_higgs_A1_A1": 4,
        "tt_201$_higgs_A1_A1_PolVec": 4,
        "tt_201$_higgs_Pi_Rho_Mixed": 10,
        "tt_201$_higgs_Pi_Pi": 4,
        "tt_201$_higgs_Pi_A1_Mixed": 4,
        "tt_201$_higgs_Pi_0A1_Mixed": 4,
        "tt_201$_higgs_A1_0A1": 4,
        "tt_201$_zttEmbed_Rho_Rho": 10,
        "tt_201$_zttEmbed_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_zttEmbed_A1_Rho": 4,
        "tt_201$_zttEmbed_A1_A1": 4,
        "tt_201$_zttEmbed_A1_A1_PolVec": 4,
        "tt_201$_zttEmbed_Pi_Rho_Mixed": 10,
        "tt_201$_zttEmbed_Pi_Pi": 4,
        "tt_201$_zttEmbed_Pi_A1_Mixed": 4,
        "tt_201$_zttEmbed_Pi_0A1_Mixed": 4,
        "tt_201$_zttEmbed_A1_0A1": 4,
        "tt_201$_jetFakes_Rho_Rho": 10,
        "tt_201$_jetFakes_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_jetFakes_A1_Rho": 4,
        "tt_201$_jetFakes_A1_A1": 4,
        "tt_201$_jetFakes_A1_A1_PolVec": 4,
        "tt_201$_jetFakes_Pi_Rho_Mixed": 10,
        "tt_201$_jetFakes_Pi_Pi": 4,
        "tt_201$_jetFakes_Pi_A1_Mixed": 4,
        "tt_201$_jetFakes_Pi_0A1_Mixed": 4,
        "tt_201$_jetFakes_A1_0A1": 4,
        "mt_ztt_201$": 1,
        "mt_fakes_201$": 1,
        "mt_murho_sig_201$": 10,
        "mt_mupi_sig_201$": 8,
        "mt_mua1_sig_201$": 4,
        "mt_mu0a1_sig_201$": 4,
        "mt_murho_ztt_201$": 10,
        "mt_mupi_ztt_201$": 8,
        "mt_mua1_ztt_201$": 4,
        "mt_mu0a1_ztt_201$": 4,
        "mt_murho_fakes_201$": 10,
        "mt_mupi_fakes_201$": 8,
        "mt_mua1_fakes_201$": 4,
        "mt_mu0a1_fakes_201$": 4,

        "et_ztt_201$": 1,
        "et_fakes_201$": 1,
        "et_murho_sig_201$": 10,
        "et_mupi_sig_201$": 8,
        "et_mua1_sig_201$": 4,
        "et_mu0a1_sig_201$": 4,
        "et_murho_ztt_201$": 10,
        "et_mupi_ztt_201$": 8,
        "et_mua1_ztt_201$": 4,
        "et_mu0a1_ztt_201$": 4,
        "et_murho_fakes_201$": 10,
        "et_mupi_fakes_201$": 8,
        "et_mua1_fakes_201$": 4,
        "et_mu0a1_fakes_201$": 4,
}

ss_bins = {

        "tt_201$_zttEmbed": 1,
        "tt_201$_jetFakes": 1,
        "tt_201$_higgs_Rho_Rho": 10,
        "tt_201$_higgs_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_higgs_A1_Rho": 4,
        "tt_201$_higgs_A1_A1": 4,
        "tt_201$_higgs_Pi_Rho_Mixed": 10,
        "tt_201$_higgs_Pi_Pi": 4,
        "tt_201$_higgs_Pi_A1_Mixed": 4,
        "tt_201$_higgs_Pi_0A1_Mixed": 4,
        "tt_201$_higgs_A1_0A1": 4,
        "tt_201$_zttEmbed_Rho_Rho": 10,
}

def MergeXBins(hist, nxbins):
  histnew = hist.Clone()
  nbins = hist.GetNbinsX()
  for i in range(1,nbins+1,nxbins):
    tot_err = ROOT.Double()
    tot = hist.IntegralAndError(i,i+nxbins-1,tot_err)
    for j in range(i,i+nxbins):
      histnew.SetBinContent(j,tot/nxbins)
      histnew.SetBinError(j,tot_err/nxbins)
  return histnew

def Symmetrise(hist,nxbins):
  histnew=hist.Clone()
  nbins = hist.GetNbinsX()
  if nbins % 2:
    print 'N X bins in 2D histogram is not even so cannot symmetrise!'
    return
  nybins = nbins/nxbins
  for i in range(1,nxbins/2+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins
      c1 = hist.GetBinContent(lo_bin_)
      c2 = hist.GetBinContent(hi_bin_)
      e1 = hist.GetBinError(lo_bin_)
      e2 = hist.GetBinError(hi_bin_)
      cnew = (c1+c2)/2
      enew = math.sqrt(e1**2 + e2**2)/2
      histnew.SetBinContent(lo_bin_,cnew)
      histnew.SetBinContent(hi_bin_,cnew)
      histnew.SetBinError(lo_bin_,enew)
      histnew.SetBinError(hi_bin_,enew)
  return histnew

def ASymmetrise(hist,hsm,hps,nxbins):
  histnew=hist.Clone()
  hsub=hsm.Clone()
  hsub.Add(hps)
  hsub.Scale(0.5)
  for i in range(1,hsub.GetNbinsX()+1): 
    histnew.SetBinContent(i,histnew.GetBinContent(i)-hsub.GetBinContent(i))
    #e = histnew.GetBinError(i) - hsm.GetBinError(i)/2 - hps.GetBinError(i)/2
    #histnew.SetBinError(i,e)
  #hsub.SetBinError(i,0.) # 0 errors as we dont want to include the subtracted components bbb uncertainties in the final uncertainty
  #histnew.Add(hsub,-1)
  nbins = hist.GetNbinsX()
  if nbins % 2:
    print 'N X bins in 2D histogram is not even so cannot symmetrise!'
    return
  nybins = nbins/nxbins
  for i in range(1,nxbins/2+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins

      mmi = hist.GetBinContent(lo_bin_)       
      mmj = hist.GetBinContent(hi_bin_)       
      smi = hsm.GetBinContent(lo_bin_) 
      smj = hsm.GetBinContent(hi_bin_)
      psi = hps.GetBinContent(lo_bin_)
      psj = hps.GetBinContent(hi_bin_)

      e_mmi = hist.GetBinError(lo_bin_)
      e_mmj = hist.GetBinError(hi_bin_)
      e_smi = hsm.GetBinError(lo_bin_)
      e_smj = hsm.GetBinError(hi_bin_)
      e_psi = hps.GetBinError(lo_bin_)
      e_psj = hps.GetBinError(hi_bin_) 

      c1_new = ( smj+psj-mmj + mmi)/2
      c2_new = ( smi+psi-mmi + mmj)/2

      e1_new = math.sqrt((e_smj+e_psj-e_mmj)**2 + e_mmi**2)/2 
      e2_new = math.sqrt((e_smi+e_psi-e_mmi)**2 + e_mmj**2)/2 

      histnew.SetBinContent(lo_bin_,c1_new)
      histnew.SetBinContent(hi_bin_,c2_new)
      histnew.SetBinError(lo_bin_,e1_new)
      histnew.SetBinError(hi_bin_,e1_new)

      #print c1_new, c2_new
      #print e1_new, e2_new    

  return histnew

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)
    year='2018'
    if '2016' in dirname: year='2016'
    if '2017' in dirname: year='2017'
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          print histo.GetName() 
          if dirname.replace(year,'201$') in cp_bins: nxbins = cp_bins[dirname.replace(year,'201$')]
          else: nxbins=1

          #for signal always symmetrise SM and PS and anti-symmetrise MM
          if 'htt125' in key.GetName() and nxbins>1:
            if '_mm_htt125' in key.GetName():
              hsm = directory.Get(key.GetName().replace('_mm_','_sm_'))
              hps = directory.Get(key.GetName().replace('_mm_','_ps_'))
              histo = ASymmetrise(histo,hsm,hps,nxbins)
            else: histo = Symmetrise(histo,nxbins)

          skip = ('data_obs' in key.GetName() or 'htt125' in key.GetName())
          rename = ('data_obs' in key.GetName())
          if rename:
            outfile.cd()
            if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
            ROOT.gDirectory.cd(dirname)
            histo.Write()
          if nxbins>1 and not skip and not rename:
            print 'rebinning for ', dirname, key.GetName()
            if '_mupi_' not in dirname and '_Pi_Pi' not in dirname and 'PolVec' not in dirname and not key.GetName().startswith('jetFakes'): histo =  MergeXBins(histo,nxbins)
            else: histo =  Symmetrise(histo,nxbins)
          if nxbins>1 and rename:
            print 'rebinning for ', dirname, key.GetName()
            if '_mupi_' not in dirname and '_Pi_Pi' not in dirname and 'PolVec' not in dirname and not key.GetName().startswith('jetFakes'): histo =  MergeXBins(histo,nxbins)
            else: histo =  Symmetrise(histo,nxbins)
            histo.SetName(histo.GetName()+'_merged')
          outfile.cd()
          if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
          ROOT.gDirectory.cd(dirname)
          print 'Writing ', dirname, histo.GetName()
          histo.Write()
          ROOT.gDirectory.cd('/')

def MergeWH(infile,outfile,dirname):
    directory = infile.Get(dirname)
    year='2018'
    if '2016' in dirname: year='2016'
    if '2017' in dirname: year='2017'
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName()).Clone()
        if 'WplusH' not in key.GetName(): continue
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          print histo.GetName()
          histo2 = directory.Get(key.GetName().replace('plus','minus'))
          print histo, histo2
          print histo.Integral(), histo2.Integral()
          histo.Add(histo2)
          print histo.Integral()
          outfile.cd()
          if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
          ROOT.gDirectory.cd(dirname)
          print 'Writing ', dirname, histo.GetName()
          histo.Write(key.GetName().replace('plus',''), ROOT.TObject.kOverwrite)
          ROOT.gDirectory.cd('/')

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge X bins')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-mergeXbins.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        #if 'murho' not in key.GetName() or 'sig' not in key.GetName(): continue
        dirname=key.GetName()
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)
        #MergeWH(output_file,output_file,dirname)

