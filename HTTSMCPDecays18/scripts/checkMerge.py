import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)

canv1 = ROOT.TCanvas()

bins = {
        "tt_201$_zttEmbed": 1,
        "tt_201$_jetFakes": 2,
        "tt_201$_higgs_Rho_Rho": 3,
        "tt_201$_higgs_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_higgs_A1_Rho": 5,
        "tt_201$_higgs_A1_A1": 6,
        "tt_201$_higgs_Pi_Rho_Mixed": 7,
        "tt_201$_higgs_Pi_Pi": 8,
        "tt_201$_higgs_Pi_A1_Mixed": 9,
        "tt_201$_higgs_Pi_0A1_Mixed": 10,
        "tt_201$_higgs_A1_0A1": 11,

        "tt_201$_zttEmbed_Rho_Rho": 3,
        "tt_201$_zttEmbed_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_zttEmbed_A1_Rho": 5,
        "tt_201$_zttEmbed_A1_A1": 6,
        "tt_201$_zttEmbed_Pi_Rho_Mixed": 7,
        "tt_201$_zttEmbed_Pi_Pi": 8,
        "tt_201$_zttEmbed_Pi_A1_Mixed": 9,
        "tt_201$_zttEmbed_Pi_0A1_Mixed": 10,
        "tt_201$_zttEmbed_A1_0A1": 11,

        "tt_201$_jetFakes_Rho_Rho": 3,
        "tt_201$_jetFakes_0A1_Rho_and_0A1_0A1": 4,
        "tt_201$_jetFakes_A1_Rho": 5,
        "tt_201$_jetFakes_A1_A1": 6,
        "tt_201$_jetFakes_Pi_Rho_Mixed": 7,
        "tt_201$_jetFakes_Pi_Pi": 8,
        "tt_201$_jetFakes_Pi_A1_Mixed": 9,
        "tt_201$_jetFakes_Pi_0A1_Mixed": 10,
        "tt_201$_jetFakes_A1_0A1": 11,
        "mt_ztt_201$": 1,
        "mt_fakes_201$": 2,
        "mt_murho_sig_201$": 3,
        "mt_mupi_sig_201$": 4,
        "mt_mua1_sig_201$": 5,
        "mt_mu0a1_sig_201$": 6,

        "mt_murho_fakes_201$": 3,
        "mt_mupi_fakes_201$": 4,
        "mt_mua1_fakes_201$": 5,
        "mt_mu0a1_fakes_201$": 6,

        "mt_murho_ztt_201$": 3,
        "mt_mupi_ztt_201$": 4,
        "mt_mua1_ztt_201$": 5,
        "mt_mu0a1_ztt_201$": 6,
}

def TwoPadSplit(split_point, gap_low, gap_high):
    upper = ROOT.TPad('upper', 'upper', 0., 0., 1., 1.)
    upper.SetBottomMargin(split_point + gap_high)
    upper.SetFillStyle(4000)
    upper.SetTicks(1)
    upper.Draw()
    lower = ROOT.TPad('lower', 'lower', 0., 0., 1., 1.)
    lower.SetTopMargin(1 - split_point + gap_low)
    lower.SetFillStyle(4000)
    lower.Draw()
    upper.cd()
    result = [upper, lower]
    return result

def checkHistogram(infile,outfile,dirname,summary,summary_fakes,summary_ztt):
    directory1 = infile.Get(dirname)
    directory2 = outfile.Get(dirname)
    count=0
    year='2018'
    if '2017' in dirname: year = '2017'
    if '2016' in dirname: year = '2016'
    for key in ['EmbedZTT','jetFakes','TTT','VVT','ZL','data_obs']:
        count+=1
        if key == 'data_obs' and 'ztt' not in dirname and 'fakes' not in dirname and 'Fakes' not in dirname: continue
        histo1 = directory1.Get(key)
        if key == 'data_obs': histo2 = directory2.Get(key+'_merged')
        else: histo2 = directory2.Get(key)
        if key == 'data_obs' and not (isinstance(histo2,ROOT.TH1D) or isinstance(histo2,ROOT.TH1F)): histo2 = directory2.Get(key)
        if histo1.Integral() <=0: continue
        ks_score = histo1.KolmogorovTest(histo2)
        print dirname, key, ks_score
        rhisto1 = histo1.Clone()
        rhisto2 = histo2.Clone()
        for i in range(1, histo1.GetNbinsX()+1):
          c1 = histo1.GetBinContent(i)
          c2 = histo2.GetBinContent(i)
          e1 = histo1.GetBinError(i)
          e2 = histo2.GetBinError(i)
          if c2>0:  
            rhisto1.SetBinContent(i,c1/c2)
            rhisto1.SetBinError(i,e1/c2)
            rhisto2.SetBinContent(i,c2/c2)
            rhisto2.SetBinError(i,e2/c2)
          else:
            rhisto1.SetBinContent(i,0)
            rhisto1.SetBinError(i,0)
            rhisto2.SetBinContent(i,0)
            rhisto2.SetBinError(i,0)
        rhisto1.SetLineColor(ROOT.kRed)
        rhisto1.SetStats(0)
        rhisto2.SetStats(0)
        histo1.SetLineColor(ROOT.kRed)
        histo1.SetStats(0)
        histo2.SetStats(0)
        histo1.GetXaxis().SetLabelSize(0)
        histo2.GetXaxis().SetLabelSize(0)
        canv1.cd()
        pads[0].cd()
        histo1.GetXaxis().SetTitle('')
        histo1.GetYaxis().SetTitle('Events')
        histo1.SetTitle('%(dirname)s %(key)s' % vars())
        histo1.Draw()
        histo2.Draw('same')
        histo1.Draw('same')
        leg = ROOT.TLegend(0.60,0.6,0.85,0.85)
        leg.AddEntry(histo1, 'un-merged', 'lep')
        leg.AddEntry(histo2, 'merged', 'lep')
        o = ROOT.TH1D()
        leg.AddEntry(o, "KS-score = %.4f" % ks_score, "")
        leg.Draw()
        pads[1].cd()
        rhisto1.GetXaxis().SetTitle('bin')
        rhisto1.GetYaxis().SetTitle('Ratio')
        rhisto1.SetTitle('')
        rhisto1.Draw()
        rhisto2.Draw('same')
        rhisto1.Draw('same')
        canv1.Print('mergeCheck_%(dirname)s_%(key)s.pdf' % vars())
        xbin = bins[dirname.replace(year,'201$')]
        if 'jetFakes' in dirname or 'fakes' in dirname: summary_fakes.SetBinContent(xbin-2,count,ks_score)
        elif 'ztt' in dirname: summary_ztt.SetBinContent(xbin-2,count,ks_score)
        else: summary.SetBinContent(xbin-2,count,ks_score)
pads=TwoPadSplit(0.29,0.01,0.01)

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge X bins')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-mergeXbins.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename) 

if '_mt.' in filename:
  summary = ROOT.TH2D('summary','',4,3,7,6,0,6)
else:
  summary = ROOT.TH2D('summary','',9,3,12,6,0,6)

summary.SetStats(0)
summary.SetTitle('')
summary.GetYaxis().SetTitle('Process')
summary.GetXaxis().SetTitle('category bin number')
summary.GetZaxis().SetTitle('KS-score')
summary.GetZaxis().SetRangeUser(0,1)

#ROOT.gStyle.SetPalette(1)
summary.GetYaxis().SetTitleOffset(1.1)
summary.GetYaxis().SetBinLabel(1,'EmbedZTT')
summary.GetYaxis().SetBinLabel(2,'jetFakes')
summary.GetYaxis().SetBinLabel(3,'TT')
summary.GetYaxis().SetBinLabel(4,'VV')
summary.GetYaxis().SetBinLabel(5,'ZL')
summary.GetYaxis().SetBinLabel(6,'data')
summary.GetXaxis().CenterLabels()
summary.GetXaxis().CenterLabels()
summary.GetXaxis().SetNdivisions(summary.GetNbinsX()+1)

summary_fakes=summary.Clone()
summary_ztt=summary.Clone()

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        #if 'other' in dirname or 'fakes' in dirname or 'ztt' in dirname or 'Fake' in dirname: continue
        if 'other' in dirname: continue

        checkHistogram(original_file,output_file,dirname,summary,summary_fakes,summary_ztt)

year='2018'
if '2016' in filename: year='2016'
if '2017' in filename: year='2017'
chan='mt'
if '_tt.' in filename: chan='tt'

canv2 = ROOT.TCanvas()
canv2.SetRightMargin(0.2)
summary.GetYaxis().SetRangeUser(0,summary.GetYaxis().GetBinLowEdge(summary.GetNbinsY()))
summary.Draw('colz')

canv2.Print('mergeCheck_summary_%(chan)s_%(year)s.pdf' % vars())

summary_fakes.Draw('colz')

canv2.Print('mergeCheck_summary_fakes_%(chan)s_%(year)s.pdf' % vars())

summary_ztt.Draw('colz')

canv2.Print('mergeCheck_summary_ztt_%(chan)s_%(year)s.pdf' % vars())
