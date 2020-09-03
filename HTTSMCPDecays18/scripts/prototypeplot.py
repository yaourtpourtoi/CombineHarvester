import numpy as np
import argparse
import ROOT
import math
import plotting2 as plot

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which binning is to be optimized')
args = parser.parse_args()
filename = args.file

f = ROOT.TFile(filename)

sig_threshold = 4.
bbb_max_threshold = 0.8
bbb_ave_threshold = 0.1

hists = []
soverb_hists = []
data_hists = []
bkg_hists = []
sm_hists = []
ps_hists = []
best_hists = []

f2 = ROOT.TFile('shapes_sm_bestfit.root')
f3 = ROOT.TFile('shapes_ps_bestfit.root')

for key in f.GetListOfKeys():
    if isinstance(f.Get(key.GetName()),ROOT.TDirectory):

      cat=key.GetName() 

      # ignore categories not sensitive to CP
      if '_1_' in cat or '_2_' in cat or 'prefit' in cat: continue
      if '_mt_' in cat: continue


      bkg = f.Get('%(cat)s/TotalBkg' % vars())
      proc = f.Get('%(cat)s/TotalProcs' % vars())

      sig = f.Get('%(cat)s/TotalSig' % vars())
      sm = f2.Get('%(cat)s/TotalSig' % vars())
      ps = f3.Get('%(cat)s/TotalSig' % vars())

      data = f.Get('%(cat)s/data_obs' % vars())

      sm_hists.append(sm.Clone())
      ps_hists.append(ps.Clone())
      best_hists.append(sig.Clone())

      x = ps.Clone()

      for i in range(1,sm.GetNbinsX()+1):
        ave = sm.GetBinContent(i)+ps.GetBinContent(i)
        diff = sm.GetBinContent(i)-ps.GetBinContent(i)
        #x.SetBinContent(i,diff/math.sqrt(ave))
        x.SetBinContent(i,diff/ave)

      hists.append(x.Clone())

      y = bkg.Clone()
      for i in range(1,sm.GetNbinsX()+1):
        ave = (sm.GetBinContent(i)+ps.GetBinContent(i))/2
        b = bkg.GetBinContent(i)
        #y.SetBinContent(i,ave/(ave+b))
        y.SetBinContent(i,ave/b)

      soverb_hists.append(y.Clone())

      data_hists.append(data.Clone())
      bkg_hists.append(bkg.Clone())



bins = np.array([-0.3,-0.1,-0.05,0.05,0.1,0.3]) # 5 bin option
bins = np.array([-0.3,-0.07,0.07,0.3]) # 3 bin option
bins = np.array([-0.3,-0.09,0.,0.09,0.3]) # 4 bin option
#bins = np.array([-0.3,-0.10,-0.06,0.06,0.10,0.3]) # 4 bin option
nbins=len(bins)-1
h = ROOT.TH1D('h','',nbins, bins)
h_sm = ROOT.TH1D('h_sm','',nbins, bins)
h_ps = ROOT.TH1D('h_ps','',nbins, bins)
h_best = ROOT.TH1D('h_best','',nbins, bins)
h_data = ROOT.TH1D('h_data','',nbins, bins)
h_bkg = ROOT.TH1D('h_bkg','',nbins, bins)

max_val=bins[-1]
min_val=bins[0]

def UpdateBinContent(hist, newc, newe, val, w):
  bini = hist.FindBin(val)
  c = hist.GetBinContent(bini)
  e = hist.GetBinError(bini)
  hist.SetBinContent(bini, c+newc*w)
  hist.SetBinError(bini, math.sqrt(e**2+(newe*w)**2))
  return hist

for x in range(0,len(hists)):
  h = hists[x]
  wts = soverb_hists[x]
  sm = sm_hists[x]
  ps = ps_hists[x]
  best = best_hists[x]
  bkg = bkg_hists[x]
  data = data_hists[x]
  for i in range(1,h.GetNbinsX()+1):
    c = h.GetBinContent(i)

    sm_newc = sm.GetBinContent(i)
    sm_newe = sm.GetBinError(i)

    ps_newc = ps.GetBinContent(i)
    ps_newe = ps.GetBinError(i)
    best_newc = best.GetBinContent(i)
    best_newe = best.GetBinError(i)


    data_newc = data.GetBinContent(i)
    data_newe = data.GetBinError(i)
    bkg_newc = bkg.GetBinContent(i)
    bkg_newe = bkg.GetBinError(i)

    w=1.
    w=wts.GetBinContent(i)

    h_sm = UpdateBinContent(h_sm, sm_newc, sm_newe, c, w)
    h_ps = UpdateBinContent(h_ps, ps_newc, ps_newe, c, w)
    h_best = UpdateBinContent(h_best, best_newc, best_newe, c, w)

    h_bkg = UpdateBinContent(h_bkg, bkg_newc, bkg_newe, c, w)
    h_data = UpdateBinContent(h_data, data_newc, data_newe, c, w)

    h.Fill(c)
    if max_val==-1 or c>max_val: max_val=c
    if min_val==-1 or c<min_val: min_val=c

c = ROOT.TCanvas()
h_ps.SetLineColor(ROOT.kRed)

def Subtract(h1,h2):
  for i in range(1,h1.GetNbinsX()+1):
    diff = h1.GetBinContent(i) - h2.GetBinContent(i)
    h1.SetBinContent(i,diff)
  return h1

#h_data.Add(h_bkg,-1)
h_data = Subtract(h_data,h_bkg)
h_bkg = Subtract(h_bkg,h_bkg)


h_data.SetLineColor(ROOT.kBlack)
h_bkg.SetLineColor(ROOT.kPink)

hs = ROOT.THStack()
hs.Add(h_sm)
hs.Add(h_ps)
hs.Add(h_data)
hs.Add(h_bkg)

#for i in range(1,h_bkg.GetNbinsX()+1):
#  h_bkg.SetBinContent(i,h_data.GetBinContent(i))

plot.propoganda_plot(h_sm,h_ps,h_best,h_bkg,h_data,'propoganda')

print h_best.GetBinContent(1), h_sm.GetBinContent(1)

