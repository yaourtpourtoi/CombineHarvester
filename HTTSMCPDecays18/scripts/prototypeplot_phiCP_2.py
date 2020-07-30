import numpy as np
import argparse
import ROOT
import math
import plotting2 as plot

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which binning is to be optimized')
args = parser.parse_args()
filename = args.file

mode=1
# 2 = rho+rho only
# 3 = pi+rho only
# 4 = mu+rho only

f = ROOT.TFile(filename)

f2 = ROOT.TFile('shapes_sm_bestfit.root')
f3 = ROOT.TFile('shapes_ps_bestfit.root')

new_sm_hist = ROOT.TH1D('new_sm_hist', '', 4,0,4)
new_ps_hist = ROOT.TH1D('new_ps_hist', '', 4,0,4)
new_bkg_hist = ROOT.TH1D('new_bkg_hist', '', 4,0,4)
new_data_hist = ROOT.TH1D('new_data_hist', '', 4,0,4)

print new_bkg_hist.GetBinError(1)

for i in range(1,new_bkg_hist.GetNbinsX()):
  new_bkg_hist.SetBinContent(i,0.)
  new_bkg_hist.SetBinError(i,0.)
  new_data_hist.SetBinContent(i,0.)
  new_data_hist.SetBinError(i,0.)
  new_sm_hist.SetBinContent(i,0.)
  new_sm_hist.SetBinError(i,0.)
  new_ps_hist.SetBinContent(i,0.)
  new_ps_hist.SetBinError(i,0.)

#def PhaseShiftMT():

for key in f.GetListOfKeys():
    if isinstance(f.Get(key.GetName()),ROOT.TDirectory):

      cat=key.GetName() 

      # ignore categories not sensitive to CP
      if ('_3_' in cat or ('_7_' in cat and '_tt_' in cat)) or 'prefit' in cat or '_1_' in cat or '_2_' in cat: continue

      print cat

      bkg = f.Get('%(cat)s/TotalBkg' % vars())
      proc = f.Get('%(cat)s/TotalProcs' % vars())

      sig = f.Get('%(cat)s/TotalSig' % vars())
      sm = f2.Get('%(cat)s/TotalSig' % vars())
      ps = f3.Get('%(cat)s/TotalSig' % vars())

      data = f.Get('%(cat)s/data_obs' % vars())

      if 'mt' in cat and '_4_' in cat:
        data.Rebin(2)
        bkg.Rebin(2)
        sm.Rebin(2)
        ps.Rebin(2)
        sig.Rebin(2)

      aves=[]
      for i in range(1,sm.GetNbinsX()+1,4):
       ave_s_p=0. 
       for j in range(1,4):
        ave = sm.GetBinContent(i)+ps.GetBinContent(i)
        diff = sm.GetBinContent(i)-ps.GetBinContent(i)
        ave_s_p+=abs(diff)/ave
       ave_s_p/=4
       aves.append(ave_s_p)
      print aves

      ave_s_p=0.
      for i in range(1,sm.GetNbinsX()+1):
        ave = sm.GetBinContent(i)+ps.GetBinContent(i)
        diff = sm.GetBinContent(i)-ps.GetBinContent(i)
        ave_s_p+=abs(diff)/ave
      ave_s_p/=sm.GetNbinsX()

      for i in range(1,sm.GetNbinsX()+1):
        bini=  int(i-math.floor(i/4)*4 )
        bin_key = int(math.floor((i-1)/4))
        if bini==0: bini=4
        mini = int(i-bini)+1
        maxi = int(mini+9)
        ave = (sm.Integral(mini,maxi)+ps.Integral(mini,maxi))/2
        #ave = sm.Integral(mini,maxi)
        b = bkg.Integral(mini,maxi)
        s_b=ave/(b+ave)
        ave_s_p = aves[bin_key]
        w = s_b*ave_s_p
        if '_mt_' in cat or ('_tt_' in cat and ('_5_' in cat or '_6_' in cat or '_9_' in cat or '_11_' in cat)):
          bini+=2
          if bini>4: bini-=4
        new_sm_hist.SetBinContent(bini, new_sm_hist.GetBinContent(bini)+w*sm.GetBinContent(i))
        new_ps_hist.SetBinContent(bini, new_ps_hist.GetBinContent(bini)+w*ps.GetBinContent(i))
        new_bkg_hist.SetBinContent(bini, new_bkg_hist.GetBinContent(bini)+w*bkg.GetBinContent(i))
        new_data_hist.SetBinContent(bini, new_data_hist.GetBinContent(bini)+w*data.GetBinContent(i))

        new_sm_hist.SetBinError(bini, math.sqrt(new_sm_hist.GetBinError(bini)**2+(w*sm.GetBinError(i))**2))
        new_ps_hist.SetBinError(bini, math.sqrt(new_ps_hist.GetBinError(bini)**2+(w*ps.GetBinError(i))**2))
        new_bkg_hist.SetBinError(bini, math.sqrt(new_bkg_hist.GetBinError(bini)**2+(w*bkg.GetBinError(i))**2))
        new_data_hist.SetBinError(bini, math.sqrt(new_data_hist.GetBinError(bini)**2+(w*data.GetBinError(i))**2))


c = ROOT.TCanvas()

def Subtract(h1,h2):
  for i in range(1,h1.GetNbinsX()+1):
    diff = h1.GetBinContent(i) - h2.GetBinContent(i)
    h1.SetBinContent(i,diff)
  return h1



#new_data_hist.Add(new_bkg_hist,-1)
new_data_hist = Subtract(new_data_hist,new_bkg_hist)
new_bkg_hist = Subtract(new_bkg_hist,new_bkg_hist)

plot_name = 'propoganda_phicp_other'

fout = ROOT.TFile('temp2.root','RECREATE')
fout.cd()
new_data_hist.Write()
new_sm_hist.Write()

plot.propoganda_plot_phicp(new_sm_hist,new_ps_hist,new_sm_hist,new_bkg_hist,new_data_hist,plot_name, 5)


