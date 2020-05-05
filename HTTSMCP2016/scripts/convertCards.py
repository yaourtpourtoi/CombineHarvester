import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

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
  hsub.Scale(0.5) # not needed if we scale MM to 2*SM cross section!!!!!!
  for i in range(1,hsub.GetNbinsX()+1): hsub.SetBinError(i,0.) # 0 errors as we dont want to include the subtracted components bbb uncertainties in the final uncertainty
  histnew.Add(hsub,-1)
  #canv = ROOT.TCanvas()
  #histnew.Draw()
  #canv.Print('test1.pdf')
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
      c1 = histnew.GetBinContent(lo_bin_)
      c2 = histnew.GetBinContent(hi_bin_)
      e1 = histnew.GetBinError(lo_bin_)
      e2 = histnew.GetBinError(hi_bin_)
      cnew = (abs(c1)+abs(c2))/2
      enew = math.sqrt(e1**2 + e2**2)/2
      c1_new = cnew*c1/abs(c1)
      c2_new = cnew*c2/abs(c2)
      histnew.SetBinContent(lo_bin_,c1_new)
      histnew.SetBinContent(hi_bin_,c2_new)
      #print lo_bin_, hi_bin_, c1, c2, c1_new, c2_new
      histnew.SetBinError(lo_bin_,enew)
      histnew.SetBinError(hi_bin_,enew)
  #histnew.Draw()
  #canv.Print('test2.pdf')
  histnew.Add(hsub)
  return histnew

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)
    year='2018'
    if '2016' in dirname: year='2016'
    if '2017' in dirname: year='2017'
    # in first loop we symmetrise all templates except for mm ones
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F): 
          nxbins=12
          skip = ('data_obs' in key.GetName())
          asymm = ('ggH_mm' in key.GetName() or 'ggHmm' in key.GetName()  or 'qqHmm' in key.GetName()  or 'qqH_mm' in key.GetName() or 'WHmm' in key.GetName() or 'WH_mm' in key.GetName()  or 'ZHmm' in key.GetName() or 'ZH_mm' in key.GetName())
          if nxbins>1 and not skip:
            print 'rebinning for ', dirname, key.GetName()
            if not asymm: histo =  Symmetrise(histo,nxbins)
          outfile.cd()
          if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
          ROOT.gDirectory.cd(dirname)
          histo.Write()
          ROOT.gDirectory.cd('/')

    # now loop to antisymmetrise mm templates
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          if not ('ggH_mm' in key.GetName() or 'ggHmm' in key.GetName()  or 'qqHmm' in key.GetName()  or 'qqH_mm' in key.GetName() or 'WHmm' in key.GetName() or 'WH_mm' in key.GetName()  or 'ZHmm' in key.GetName() or 'ZH_mm' in key.GetName()): continue
          nxbins=12
          if nxbins>1 and not skip:
            print 'rebinning for ', dirname, key.GetName()
            sm_name = key.GetName()
            ps_name = key.GetName()
            to_asymm = ['ggHmm', 'ggH_mm', 'qqHmm','qqH_mm','WHmm','WH_mm','ZHmm','ZH_mm']
            for x in to_asymm:
              if x in key.GetName():
                sm_name = sm_name.replace('mm','sm')
                ps_name = ps_name.replace('mm','ps')
                break
            hsm = directory.Get(sm_name)
            hps = directory.Get(ps_name)
            histo =  ASymmetrise(histo,hsm,hps,nxbins)
          outfile.cd()
          if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
          ROOT.gDirectory.cd(dirname)
          histo.Write()
          ROOT.gDirectory.cd('/')


parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge X bins')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-symm.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        #if not ('boosted' in dirname and 'tightmjj' in dirname): continue
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)

