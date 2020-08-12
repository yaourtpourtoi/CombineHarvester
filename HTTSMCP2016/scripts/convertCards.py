import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)

#canv = ROOT.TCanvas()
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

def WeightedAve(h1,h2,h3):
  hout = h1.Clone()
  for i in range(1,hout.GetNbinsX()+1):

    if h1.GetBinError(i) > 0: w1 = 1./h1.GetBinError(i)**2
    else: w1 = 0.
    if h2.GetBinError(i) > 0: w2 = 1./h2.GetBinError(i)**2
    else: w2 = 0.
    if h3.GetBinError(i) > 0: w3 = 1./h3.GetBinError(i)**2
    else: w3 = 0.


    c1 = h1.GetBinContent(i)
    c2 = h2.GetBinContent(i)
    c3 = h3.GetBinContent(i)

    if (w1+w2+w3) != 0: 
      c = (w1*c1+w2*c2+w3*c3)/(w1+w2+w3)
      u = 1./math.sqrt(w1+w2+w3)
    else:
      c = 0.
      u = 0.

    hout.SetBinContent(i,c)
    hout.SetBinError(i,u)

  return hout

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

def ASymmetrise(hist,hsm,hps,nxbins,vbf=False):
  histnew=hist.Clone()
  hsub=hsm.Clone()
  hsub.Add(hps)
  if vbf: hsub.Scale(0.5)
  for i in range(1,hsub.GetNbinsX()+1): hsub.SetBinError(i,0.) # 0 errors as we dont want to include the subtracted components bbb uncertainties in the final uncertainty
  histnew.Add(hsub,-1)
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

      if c1*c2 < 0:
        c1_new = cnew*c1/abs(c1)
        c2_new = cnew*c2/abs(c2)
      else:
        if c1 < c2:
          c1_new = -cnew
          c2_new = cnew
        else: 
          c1_new = cnew
          c2_new = -cnew
        
      histnew.SetBinContent(lo_bin_,c1_new)
      histnew.SetBinContent(hi_bin_,c2_new)
      histnew.SetBinError(lo_bin_,enew)
      histnew.SetBinError(hi_bin_,enew)
  histnew.Add(hsub)
  return histnew


def makeVBFUncert(h, h_ph, h_jhu):
  nxbins=12
  nybins = h.GetNbinsX()/nxbins

  h_merge_ph = ROOT.TH1D('h_merge_ph','',nxbins,0,nxbins)
  h_merge_jhu = ROOT.TH1D('h_merge_jhu','',nxbins,0,nxbins)
  for i in range(1,nxbins+1):
    c1=0.
    c2=0.
    e1=0.
    e2=0.
    for j in range(0,nybins):
      c1+= h_ph.GetBinContent(i+j*nxbins) 
      c2+= h_jhu.GetBinContent(i+j*nxbins) 
      #e1+= h_ph.GetBinError(i+j*nxbins)**2
      #e2+= h_jhu.GetBinError(i+j*nxbins)**2
    h_merge_ph.SetBinContent(i,c1)
    h_merge_jhu.SetBinContent(i,c2)
    h_merge_ph.SetBinError(i,math.sqrt(e1))    
    h_merge_jhu.SetBinError(i,math.sqrt(e2))    
  

  h_merge_ph.Scale(1/h_merge_ph.Integral())
  h_merge_jhu.Scale(1/h_merge_jhu.Integral())
  h_merge_ph.Divide(h_merge_jhu)

  #h_merge_ph.SetStats(0)
  #h_merge_ph.GetXaxis().SetTitle('bin')
  #h_merge_ph.GetYaxis().SetTitle('NLO / LO')
  #h_merge_ph.Draw()
  #canv.Print('test.pdf')
  #exit()

  hnew_up = h.Clone()
  hnew_down = h.Clone()
  hnew_up.SetName(h.GetName()+'_VBF_NLOvsLOUp')
  hnew_down.SetName(h.GetName()+'_VBF_NLOvsLODown')

  for i in range(1,nxbins+1):
    up=h_merge_ph.GetBinContent(i)
    down=1./up
    for j in range(0,nybins):
      c_up  = hnew_up.GetBinContent(i+j*nxbins)*up
      c_down  = hnew_down.GetBinContent(i+j*nxbins)*down
      e_up = hnew_up.GetBinError(i+j*nxbins)*up
      e_down = hnew_down.GetBinError(i+j*nxbins)*down
      hnew_up.SetBinContent(i+j*nxbins,c_up)
      hnew_down.SetBinContent(i+j*nxbins,c_down)
      hnew_up.SetBinError(i+j*nxbins,math.sqrt(e_up))
      hnew_down.SetBinError(i+j*nxbins,math.sqrt(e_down))

  return (hnew_up, hnew_down)   

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname, reweight=False):
    directory = infile.Get(dirname)
    year='2018'
    if '2016' in dirname: year='2016'
    if '2017' in dirname: year='2017'

    vbf_jhu = directory.Get('qqH_sm_htt125')
    zh_jhu = directory.Get('ZH_sm_htt125')
    wh_jhu = directory.Get('WH_sm_htt125')

    vbf_ph = directory.Get('qqH_htt125')
    zh_ph = directory.Get('ZH_htt125')
    wh_ph = directory.Get('WH_htt125')

    vbf_sf=1.
    zh_sf=1.
    wh_sf=1.

    if(vbf_jhu.Integral()>0): vbf_sf = vbf_ph.Integral()/vbf_jhu.Integral()
    if(zh_jhu.Integral()>0): zh_sf = zh_ph.Integral()/zh_jhu.Integral()
    if(wh_jhu.Integral()>0): wh_sf = wh_ph.Integral()/wh_jhu.Integral()


    # if we are combining reweighted templates with standard ones then we perform this step first
    if reweight:
      for key in directory.GetListOfKeys():
          histo = directory.Get(key.GetName())
          if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):

            cps = ['sm','ps','mm']
            for i in cps:
              if i == 'sm':
                o1 = 'ps'
                o2 = 'mm'
              if i == 'ps':
                o1 = 'sm'
                o2 = 'mm'
              if i == 'mm':
                o1 = 'sm'
                o2 = 'ps'
              if 'ggH_%s_htt' % i in key.GetName() and 'reweightedto' not in key.GetName():
                h2 = directory.Get(key.GetName().replace(i+'_htt125',o1+'_htt125_reweightedto_'+i))
                h3 = directory.Get(key.GetName().replace(i+'_htt125',o2+'_htt125_reweightedto_'+i))
                histo=WeightedAve(histo,h2,h3)
                break
            outfile.cd()
            if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
            ROOT.gDirectory.cd(dirname)
            histo.Write("",ROOT.TObject.kOverwrite)
            ROOT.gDirectory.cd('/')
      
      directory = outfile.Get(dirname)

    to_write = []

    print(directory)
    # in second loop we symmetrise all templates except for mm ones
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          if 'qqH_sm_' in key.GetName() or 'qqH_mm_' in key.GetName() or 'qqH_ps_' in key.GetName():
            print 'applying SF = ', vbf_sf, ' to ', key.GetName()
            histo.Scale(vbf_sf)
          if 'WH_sm_' in key.GetName() or 'WH_mm_' in key.GetName() or 'WH_ps_' in key.GetName():
            print 'applying SF = ', wh_sf, ' to ', key.GetName()
            histo.Scale(wh_sf)
          if 'ZH_sm_' in key.GetName() or 'ZH_mm_' in key.GetName() or 'ZH_ps_' in key.GetName():
            print 'applying SF = ', zh_sf, ' to ', key.GetName()
            histo.Scale(zh_sf)
          nxbins=12
          skip = ('data_obs' in key.GetName()) or not ('dijet' in dirname or 'vbf' in dirname)
          asymm = (('ggH_mm' in key.GetName() or 'ggHmm' in key.GetName()  or 'qqHmm' in key.GetName()  or 'qqH_mm' in key.GetName() or 'WHmm' in key.GetName() or 'WH_mm' in key.GetName()  or 'ZHmm' in key.GetName() or 'ZH_mm' in key.GetName() or 'reweightedto_mm' in key.GetName())) and 'reweightedto_sm' not in key.GetName() and 'reweightedto_ps' not in key.GetName()
          if nxbins>1 and not skip and not asymm:
            print 'rebinning for ', dirname, key.GetName()
            if not asymm: histo =  Symmetrise(histo,nxbins)
          to_write.append([dirname,histo])

    for x in to_write:
      outfile.cd()
      if not ROOT.gDirectory.GetDirectory(x[0]): ROOT.gDirectory.mkdir(x[0])
      ROOT.gDirectory.cd(x[0])
      x[1].Write("",ROOT.TObject.kOverwrite)
      ROOT.gDirectory.cd('/')

    directory = outfile.Get(dirname)
    # now loop to antisymmetrise mm templates
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          nxbins=12
          skip = ('data_obs' in key.GetName()) or not ('dijet' in dirname or 'vbf' in dirname)
          if not ((('ggH_mm' in key.GetName() or 'ggHmm' in key.GetName()  or 'qqHmm' in key.GetName()  or 'qqH_mm' in key.GetName() or 'WHmm' in key.GetName() or 'WH_mm' in key.GetName()  or 'ZHmm' in key.GetName() or 'ZH_mm' in key.GetName() or 'reweightedto_mm' in key.GetName())) and 'reweightedto_sm' not in key.GetName() and 'reweightedto_ps' not in key.GetName()): continue
          nxbins=12
          if nxbins>1 and not skip:
            print 'asymm rebinning for ', dirname, key.GetName()
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
            isVBF=False
            if 'qqH' in key.GetName() or 'WH' in key.GetName() or 'ZH' in key.GetName(): isVBF=True
            histo =  ASymmetrise(histo,hsm,hps,nxbins,isVBF)
          to_write.append([dirname,histo])
          # after all symmetrisations are done add the LO vs NLO uncertainty
          if ('dijet' in dirname or 'vbf' in dirname) and (key.GetName() == 'qqH_sm_htt125' or key.GetName() == 'qqH_mm_htt125' or key.GetName() == 'qqH_ps_htt125'):
            vbf_jhu = directory.Get('qqH_sm_htt125').Clone()
            vbf_ph = directory.Get('qqH_htt125').Clone() 
            (h_up, h_down) = makeVBFUncert(histo, vbf_ph, vbf_jhu)
            to_write.append([dirname,h_up])
            to_write.append([dirname,h_down])
            print 'making uncert for', dirname, key.GetName()

    for x in to_write:
      outfile.cd()
      if not ROOT.gDirectory.GetDirectory(x[0]): ROOT.gDirectory.mkdir(x[0])
      ROOT.gDirectory.cd(x[0])
      x[1].Write("",ROOT.TObject.kOverwrite)
      ROOT.gDirectory.cd('/')

    to_write = []

    directory = outfile.Get(dirname)
    # now loop to make VBF NLO vs LO uncertainty templates
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
          # after all symmetrisations are done add the LO vs NLO uncertainty
          if ('dijet' in dirname or 'vbf' in dirname) and (key.GetName() == 'qqH_sm_htt125' or key.GetName() == 'qqH_mm_htt125' or key.GetName() == 'qqH_ps_htt125'):
            vbf_jhu = directory.Get('qqH_sm_htt125').Clone()
            vbf_ph = directory.Get('qqH_htt125').Clone()
            (h_up, h_down) = makeVBFUncert(histo, vbf_ph, vbf_jhu)
            to_write.append([dirname,h_up])
            to_write.append([dirname,h_down])
            print 'making uncert for', dirname, key.GetName()

    for x in to_write:
      outfile.cd()
      if not ROOT.gDirectory.GetDirectory(x[0]): ROOT.gDirectory.mkdir(x[0])
      ROOT.gDirectory.cd(x[0])
      x[1].Write("",ROOT.TObject.kOverwrite)
      ROOT.gDirectory.cd('/')

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge X bins')
parser.add_argument('--reweight', '-r', action='store_true', help= 'Combine reweighted templates with standard ones')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-symm.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        #if not ('dijet' in dirname or 'vbf' in dirname): continue
        print 'performing (anti-)symmetrisation for cat:', dirname
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname,reweight=args.reweight)

