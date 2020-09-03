import numpy as np
import argparse
import ROOT
import math
import plotting2 as plot

parser = argparse.ArgumentParser()
parser.add_argument('--mode', '-m', default=1, type=int  ,help= 'File from which binning is to be optimized')
args = parser.parse_args()

mode=args.mode
# 1 = rho+rho + pi+rho + mu+pi
# 2 = rho+rho only
# 3 = pi+rho only
# 4 = mu+rho only
# 5 = others

if mode ==1:
  f1 = ROOT.TFile("shapes_prop_best_sm.root")
  f2 = ROOT.TFile('shapes_prop_best_ps.root')
if mode ==2:
  f1 = ROOT.TFile("shapes_prop_tt_3_sm.root")
  f2 = ROOT.TFile('shapes_prop_tt_3_ps.root')
if mode ==3:
  f1 = ROOT.TFile("shapes_prop_tt_7_sm.root")
  f2 = ROOT.TFile('shapes_prop_tt_7_ps.root')
if mode ==4:
  f1 = ROOT.TFile("shapes_prop_mt_3_sm.root")
  f2 = ROOT.TFile('shapes_prop_mt_3_ps.root')
if mode ==5:
  f1 = ROOT.TFile("shapes_prop_worst_sm.root")
  f2 = ROOT.TFile('shapes_prop_worst_ps.root')


sm_hist = f1.Get('postfit/TotalSig')
ps_hist = f2.Get('postfit/TotalSig')
bkg_hist = f1.Get('postfit/TotalBkg')
data_hist = f1.Get('postfit/data_obs')


def Subtract(h1,h2):
  for i in range(1,h1.GetNbinsX()+1):
    diff = h1.GetBinContent(i) - h2.GetBinContent(i)
    h1.SetBinContent(i,diff)
  return h1

data_hist = Subtract(data_hist,bkg_hist)
bkg_hist = Subtract(bkg_hist,bkg_hist)

plot_name = 'phiCPWeighted'
if mode == 2: plot_name+='_rhorho'
if mode == 3: plot_name+='_pirho'
if mode == 4: plot_name+='_murho'
if mode == 5: plot_name+='_others'

fout = ROOT.TFile('temp.root','RECREATE')
fout.cd()
data_hist.Write()
sm_hist.Write()

plot.propoganda_plot_phicp(sm_hist,ps_hist,sm_hist,bkg_hist,data_hist,plot_name, mode)


