import ROOT as R
import math
COL_STORE = []

def SetTDRStyle():
    """Sets the PubComm recommended style

    Just a copy of <http://ghm.web.cern.ch/ghm/plots/MacroExample/tdrstyle.C>
    @sa ModTDRStyle() to use this style with some additional customisation.
    """
    # For the canvas:
    R.gStyle.SetCanvasBorderMode(0)
    R.gStyle.SetCanvasColor(R.kWhite)
    R.gStyle.SetCanvasDefH(600)  # Height of canvas
    R.gStyle.SetCanvasDefW(600)  # Width of canvas
    R.gStyle.SetCanvasDefX(0)    # POsition on screen
    R.gStyle.SetCanvasDefY(0)

    # For the Pad:
    R.gStyle.SetPadBorderMode(0)
    # R.gStyle.SetPadBorderSize(Width_t size = 1)
    R.gStyle.SetPadColor(R.kWhite)
    R.gStyle.SetPadGridX(False)
    R.gStyle.SetPadGridY(False)
    R.gStyle.SetGridColor(0)
    R.gStyle.SetGridStyle(3)
    R.gStyle.SetGridWidth(1)

    # For the frame:
    R.gStyle.SetFrameBorderMode(0)
    R.gStyle.SetFrameBorderSize(1)
    R.gStyle.SetFrameFillColor(0)
    R.gStyle.SetFrameFillStyle(0)
    R.gStyle.SetFrameLineColor(1)
    R.gStyle.SetFrameLineStyle(1)
    R.gStyle.SetFrameLineWidth(1)

    # For the histo:
    # R.gStyle.SetHistFillColor(1)
    # R.gStyle.SetHistFillStyle(0)
    R.gStyle.SetHistLineColor(1)
    R.gStyle.SetHistLineStyle(0)
    R.gStyle.SetHistLineWidth(1)
    # R.gStyle.SetLegoInnerR(Float_t rad = 0.5)
    # R.gStyle.SetNumberContours(Int_t number = 20)

    R.gStyle.SetEndErrorSize(2)
    # R.gStyle.SetErrorMarker(20)
    # R.gStyle.SetErrorX(0.)

    R.gStyle.SetMarkerStyle(20)

    # For the fit/function:
    R.gStyle.SetOptFit(1)
    R.gStyle.SetFitFormat('5.4g')
    R.gStyle.SetFuncColor(2)
    R.gStyle.SetFuncStyle(1)
    R.gStyle.SetFuncWidth(1)

    # For the date:
    R.gStyle.SetOptDate(0)
    # R.gStyle.SetDateX(Float_t x = 0.01)
    # R.gStyle.SetDateY(Float_t y = 0.01)

    # For the statistics box:
    R.gStyle.SetOptFile(0)
    R.gStyle.SetOptStat(0)
    # To display the mean and RMS:   SetOptStat('mr')
    R.gStyle.SetStatColor(R.kWhite)
    R.gStyle.SetStatFont(42)
    R.gStyle.SetStatFontSize(0.025)
    R.gStyle.SetStatTextColor(1)
    R.gStyle.SetStatFormat('6.4g')
    R.gStyle.SetStatBorderSize(1)
    R.gStyle.SetStatH(0.1)
    R.gStyle.SetStatW(0.15)
    # R.gStyle.SetStatStyle(Style_t style = 1001)
    # R.gStyle.SetStatX(Float_t x = 0)
    # R.gStyle.SetStatY(Float_t y = 0)

    # Margins:
    R.gStyle.SetPadTopMargin(0.05)
    R.gStyle.SetPadBottomMargin(0.13)
    R.gStyle.SetPadLeftMargin(0.16)
    R.gStyle.SetPadRightMargin(0.02)

    # For the Global title:
    R.gStyle.SetOptTitle(0)
    R.gStyle.SetTitleFont(42)
    R.gStyle.SetTitleColor(1)
    R.gStyle.SetTitleTextColor(1)
    R.gStyle.SetTitleFillColor(10)
    R.gStyle.SetTitleFontSize(0.05)
    # R.gStyle.SetTitleH(0); # Set the height of the title box
    # R.gStyle.SetTitleW(0); # Set the width of the title box
    # R.gStyle.SetTitleX(0); # Set the position of the title box
    # R.gStyle.SetTitleY(0.985); # Set the position of the title box
    # R.gStyle.SetTitleStyle(Style_t style = 1001)
    # R.gStyle.SetTitleBorderSize(2)

    # For the axis titles:
    R.gStyle.SetTitleColor(1, 'XYZ')
    R.gStyle.SetTitleFont(42, 'XYZ')
    R.gStyle.SetTitleSize(0.06, 'XYZ')
    # Another way to set the size?
    # R.gStyle.SetTitleXSize(Float_t size = 0.02)
    # R.gStyle.SetTitleYSize(Float_t size = 0.02)
    R.gStyle.SetTitleXOffset(0.9)
    R.gStyle.SetTitleYOffset(1.25)
    # R.gStyle.SetTitleOffset(1.1, 'Y'); # Another way to set the Offset

    # For the axis labels:

    R.gStyle.SetLabelColor(1, 'XYZ')
    R.gStyle.SetLabelFont(42, 'XYZ')
    R.gStyle.SetLabelOffset(0.007, 'XYZ')
    R.gStyle.SetLabelSize(0.05, 'XYZ')

    # For the axis:

    R.gStyle.SetAxisColor(1, 'XYZ')
    R.gStyle.SetStripDecimals(True)
    R.gStyle.SetTickLength(0.03, 'XYZ')
    R.gStyle.SetNdivisions(510, 'XYZ')
    R.gStyle.SetPadTickX(1)
    R.gStyle.SetPadTickY(1)

    # Change for log plots:
    R.gStyle.SetOptLogx(0)
    R.gStyle.SetOptLogy(0)
    R.gStyle.SetOptLogz(0)

    # Postscript options:
    R.gStyle.SetPaperSize(20., 20.)
    # R.gStyle.SetLineScalePS(Float_t scale = 3)
    # R.gStyle.SetLineStyleString(Int_t i, const char* text)
    # R.gStyle.SetHeaderPS(const char* header)
    # R.gStyle.SetTitlePS(const char* pstitle)

    # R.gStyle.SetBarOffset(Float_t baroff = 0.5)
    # R.gStyle.SetBarWidth(Float_t barwidth = 0.5)
    # R.gStyle.SetPaintTextFormat(const char* format = 'g')
    # R.gStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
    # R.gStyle.SetTimeOffset(Double_t toffset)
    # R.gStyle.SetHistMinimumZero(kTRUE)

    R.gStyle.SetHatchesLineWidth(5)
    R.gStyle.SetHatchesSpacing(0.05)


def ModTDRStyle(width=600, height=600, t=0.06, b=0.12, l=0.16, r=0.04):
    """Modified version of the tdrStyle

    Args:
        width (int): Canvas width in pixels
        height (int): Canvas height in pixels
        t (float): Pad top margin [0-1]
        b (float): Pad bottom margin [0-1]
        l (float): Pad left margin [0-1]
        r (float): Pad right margin [0-1]
    """
    SetTDRStyle()

    # Set the default canvas width and height in pixels
    R.gStyle.SetCanvasDefW(width)
    R.gStyle.SetCanvasDefH(height)

    # Set the default margins. These are given as fractions of the pad height
    # for `Top` and `Bottom` and the pad width for `Left` and `Right`. But we
    # want to specify all of these as fractions of the shortest length.
    def_w = float(R.gStyle.GetCanvasDefW())
    def_h = float(R.gStyle.GetCanvasDefH())

    scale_h = (def_w / def_h) if (def_h > def_w) else 1.
    scale_w = (def_h / def_w) if (def_w > def_h) else 1.

    def_min = def_h if (def_h < def_w) else def_w

    R.gStyle.SetPadTopMargin(t * scale_h)
    # default 0.05
    R.gStyle.SetPadBottomMargin(b * scale_h)
    # default 0.13
    R.gStyle.SetPadLeftMargin(l * scale_w)
    # default 0.16
    R.gStyle.SetPadRightMargin(r * scale_w)
    # default 0.02
    # But note the new CMS style sets these:
    # 0.08, 0.12, 0.12, 0.04

    # Set number of axis tick divisions
    R.gStyle.SetNdivisions(510, 'XYZ')  # default 510

    # Some marker properties not set in the default tdr style
    R.gStyle.SetMarkerColor(R.kBlack)
    R.gStyle.SetMarkerSize(1.0)

    R.gStyle.SetLabelOffset(0.007, 'YZ')
    # This is an adhoc adjustment to scale the x-axis label
    # offset when we stretch plot vertically
    # Will also need to increase if first x-axis label has more than one digit
    R.gStyle.SetLabelOffset(0.005 * (3. - 2. / scale_h), 'X')

    # In this next part we do a slightly involved calculation to set the axis
    # title offsets, depending on the values of the TPad dimensions and
    # margins. This is to try and ensure that regardless of how these pad
    # values are set, the axis titles will be located towards the edges of the
    # canvas and not get pushed off the edge - which can often happen if a
    # fixed value is used.
    title_size = 0.05
    title_px = title_size * def_min
    label_size = 0.04
    R.gStyle.SetTitleSize(title_size, 'XYZ')
    R.gStyle.SetLabelSize(label_size, 'XYZ')

    R.gStyle.SetTitleXOffset(0.5 * scale_h *
                             (1.2 * (def_h * b * scale_h - 0.6 * title_px)) /
                             title_px)
    R.gStyle.SetTitleYOffset(0.5 * scale_w *
                             (1.2 * (def_w * l * scale_w - 0.6 * title_px)) /
                             title_px)

    # Only draw ticks where we have an axis
    R.gStyle.SetPadTickX(0)
    R.gStyle.SetPadTickY(0)
    R.gStyle.SetTickLength(0.02, 'XYZ')

    R.gStyle.SetLegendBorderSize(0)
    R.gStyle.SetLegendFont(42)
    R.gStyle.SetLegendFillColor(0)
    R.gStyle.SetFillColor(0)

    R.gROOT.ForceStyle()

def OnePad():
    pad = R.TPad('pad', 'pad', 0., 0., 1., 1.)
    pad.SetTicks(1)
    pad.Draw()
    pad.cd()
    result = [pad]
    return result

def DrawCMSLogo(pad, cmsText, extraText, iPosX, relPosX, relPosY, relExtraDY, relExtraDX, extraText2='', cmsTextSize=0.8):
    """Blah
    
    Args:
        pad (TYPE): Description
        cmsText (TYPE): Description
        extraText (TYPE): Description
        iPosX (TYPE): Description
        relPosX (TYPE): Description
        relPosY (TYPE): Description
        relExtraDY (TYPE): Description
        extraText2 (str): Description
        cmsTextSize (float): Description
    
    Returns:
        TYPE: Description
    """
    pad.cd()
    cmsTextFont = 62  # default is helvetic-bold

    writeExtraText = len(extraText) > 0
    writeExtraText2 = len(extraText2) > 0
    extraTextFont = 52

    # text sizes and text offsets with respect to the top frame
    # in unit of the top margin size
    lumiTextOffset = 0.2
    # cmsTextSize = 0.8
    # float cmsTextOffset    = 0.1;  // only used in outOfFrame version

    # ratio of 'CMS' and extra text size
    extraOverCmsTextSize = 0.76

    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True

    alignY_ = 3
    alignX_ = 2
    if (iPosX / 10 == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignY_ = 1
    if (iPosX / 10 == 1):
        alignX_ = 1
    if (iPosX / 10 == 2):
        alignX_ = 2
    if (iPosX / 10 == 3):
        alignX_ = 3
    # if (iPosX == 0): relPosX = 0.14
    align_ = 10 * alignX_ + alignY_

    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()

    latex = R.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(R.kBlack)

    extraTextSize = extraOverCmsTextSize * cmsTextSize
    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if (pad_ratio < 1.):
        pad_ratio = 1.

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)

    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif (iPosX % 10 == 2):
        posX_ = l + 0.5 * (1 - l - r)
    elif (iPosX % 10 == 3):
        posX_ = 1 - r - relPosX * (1 - l - r)

    posY_ = 1 - t - relPosY * (1 - t - b)
    if not outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, cmsText)
        if writeExtraText:
            latex.SetTextFont(extraTextFont)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extraTextSize * t * pad_ratio)
            latex.DrawLatex(
                posX_+ relExtraDX * cmsTextSize * t, posY_ - relExtraDY * cmsTextSize * t, extraText)
            if writeExtraText2:
                latex.DrawLatex(
                    posX_, posY_ - 1.8 * relExtraDY * cmsTextSize * t, extraText2)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r)
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)


def PositionedLegend(width, height, pos, offset,yoffset):
    o = offset
    w = width
    h = height
    l = R.gPad.GetLeftMargin()
    t = R.gPad.GetTopMargin()
    b = R.gPad.GetBottomMargin()
    r = R.gPad.GetRightMargin()
    if pos == 1:
        return R.TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, '', 'NBNDC')
    if pos == 2:
        c = l + 0.5 * (1 - l - r)
        return R.TLegend(c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, '', 'NBNDC')
    if pos == 3:
        return R.TLegend(1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, '', 'NBNDC')
    if pos == 4:
        return R.TLegend(l + o, b + o, l + o + w, b + o + h, '', 'NBNDC')
    if pos == 5:
        c = l + 0.5 * (1 - l - r)
        return R.TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, '', 'NBNDC')
    if pos == 6:
        return R.TLegend(1 - r - o - w, b + o+yoffset, 1 - r - o, b + o+yoffset + h, '', 'NBNDC')
    if pos == 7:
        return R.TLegend(1 - o - w, 1 - t - o - h, 1 - o, 1 - t - o, '', 'NBNDC')

def DrawTitle(pad, text, align, scale=1):
    pad_backup = R.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if pad_ratio < 1.:
        pad_ratio = 1.

    textSize = 0.8
    textOffset = 0.2

    latex = R.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(R.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t
    if align == 1:
        latex.SetTextAlign(11)
    if align == 1:
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
    if align == 2:
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
    if align == 3:
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

def CreateTransparentColor(color, alpha):
    adapt = R.gROOT.GetColor(color)
    new_idx = R.gROOT.GetListOfColors().GetLast() + 1
    trans = R.TColor(
        new_idx, adapt.GetRed(), adapt.GetGreen(), adapt.GetBlue(), '', alpha)
    COL_STORE.append(trans)
    trans.SetName('userColor%i' % new_idx)
    return new_idx

def propoganda_plot(sm,ps,best,bkg,data,plot_name):

    bins = []
    for i in range(1,data.GetNbinsX()+2):
      bins.append(data.GetBinLowEdge(i))


    def ConvertToEqualBins(h):
      hnew = R.TH1D(h.GetName()+'_new','',h.GetNbinsX(),0,h.GetNbinsX())
      for i in range(1,h.GetNbinsX()+1):
        hnew.SetBinContent(i,h.GetBinContent(i))
        hnew.SetBinError(i,h.GetBinError(i))
      return hnew

    sm = ConvertToEqualBins(sm)
    ps = ConvertToEqualBins(ps)
    best = ConvertToEqualBins(best)
    bkg = ConvertToEqualBins(bkg)
    data = ConvertToEqualBins(data)
  
    data.GetXaxis().SetLabelSize(0.07) 
    data.GetYaxis().SetLabelSize(0.05) 
    data.GetXaxis().SetTitle('(S #minus P)/(S + P)')
    data.GetXaxis().SetTitleOffset(1.1)
    data.GetXaxis().SetTitleSize(0.05)
    data.GetYaxis().SetTitle('(S+P)/2B Weighted Events / bin')
    data.GetYaxis().SetTitleSize(0.05)
    for i in range(1,data.GetNbinsX()+1):
      label = ''
      print i, bins
      if i == 1: label = '< %.2f' % bins[1]
      elif i == len(bins)-1: label = '#geq %.2f' % bins[i-1]
      else: label = '%.2f #minus  %.2f' %(bins[i-1], bins[i])
      data.GetXaxis().SetBinLabel(i,label) 

    c1 = R.TCanvas()

    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    pads=OnePad()
    pads[0].cd()

    hs = R.THStack("hs","")

    data.SetMarkerStyle(20)
    data.SetLineColor(1)
    miny=0.
    maxe=0.
    for i in range(1,bkg.GetNbinsX()+1):
     e = bkg.GetBinError(i)
     if e> maxe: maxe=e
    miny=-maxe*1.4
    data.SetMinimum(miny)
    data.Draw("E")

    sm.SetLineWidth(2)
    sm.SetLineColor(R.kBlue)
    sm.SetMarkerSize(0)
    sm.SetFillStyle(0)

    ps.SetLineWidth(2)
    ps.SetLineColor(R.kRed)
    ps.SetMarkerSize(0)
    ps.SetFillStyle(0)

    best.SetLineWidth(2)
    best.SetLineStyle(2)
    best.SetLineColor(R.kBlack)
    best.SetMarkerSize(0)
    best.SetFillStyle(0)

    hs.Add(ps)
    hs.Add(sm)

    hs.Draw("nostack hist same")

    bkg.SetFillColor(CreateTransparentColor(12,0.4))
    bkg.SetLineColor(CreateTransparentColor(12,0.4))
    bkg.SetMarkerSize(0)
    bkg.SetMarkerColor(CreateTransparentColor(12,0.4))

    bkg.Draw("e2same")
    data.Draw("E same")


    DrawCMSLogo(pads[0], 'CMS', 'Preliminary', 11, 0.001, -0.07, 0.2, 1.5, '', 1.0)
    DrawTitle(pads[0], '137 fb^{-1} (13 TeV)', 3)

    #Setup legend
    legend = PositionedLegend(0.25,0.3,6,0.02,0.08)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    legend.AddEntry(data,'Data #minus Bkg.',"lep")
    legend.AddEntry(bkg,'Bkg. uncert.',"f")
    legend.AddEntry(sm,'#phi_{#tau#tau} = 0^{#circ}',"l")
    legend.AddEntry(ps,'#phi_{#tau#tau} = 90^{#circ}',"l")
    legend.Draw("same")

    c1.SaveAs(plot_name+'.pdf')

def propoganda_plot_phicp(sm,ps,best,bkg,data,plot_name,mode=1):

    title='#rho#rho + #pi#rho + #mu#rho'
    if mode == 2:
      title='#rho#rho'
    if mode == 3:
      title='#pi#rho'
    if mode == 4:
      title='#mu#rho' 
    if mode == 5:
      title = 'others'

    latex = R.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(R.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.06)


    def ConvertToEqualBins(h):
      hnew = R.TH1D(h.GetName()+'_new','',h.GetNbinsX(),0,360)
      for i in range(1,h.GetNbinsX()+1):
        hnew.SetBinContent(i,h.GetBinContent(i))
        hnew.SetBinError(i,h.GetBinError(i))
      return hnew
                                   
    #sm = ConvertToEqualBins(sm)
    #ps = ConvertToEqualBins(ps)
    #best = ConvertToEqualBins(best)    
    #bkg = ConvertToEqualBins(bkg)
    #data = ConvertToEqualBins(data)

    data.GetXaxis().SetTitleOffset(1.0)
    data.GetXaxis().SetTitleSize(0.05)
    data.GetYaxis().SetTitle('A#times S/(S+B) Weighted Events / bin')
    data.GetXaxis().SetTitle('#phi_{CP} (degrees)')
    data.GetYaxis().SetTitleSize(0.05)


    c1 = R.TCanvas()

    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    pads=OnePad()
    pads[0].cd()

    hs = R.THStack("hs","")

    data.SetMarkerStyle(20)
    data.SetLineColor(1)
    miny=0.
    maxe=0.
    for i in range(1,bkg.GetNbinsX()+1):
     e = bkg.GetBinError(i)
     if e> maxe: maxe=e
    miny=-maxe*1.4
    for i in range(1,data.GetNbinsX()+1): 
      x=data.GetBinContent(i) - data.GetBinError(i)
      if x < miny: miny=x 
    if miny<data.GetMinimum(): data.SetMinimum(miny)
    if mode ==5: data.SetMaximum(data.GetMaximum()*1.8)
    data.Draw("E")

    col_sm = R.TColor().GetColor('#253494') 
    col_ps = R.TColor().GetColor('#006837') 
    col_mix = R.TColor().GetColor('#D62839') 

    sm.SetLineWidth(3)
    sm.SetLineColor(col_sm)
    sm.SetMarkerSize(0)
    sm.SetFillStyle(0)

    ps.SetLineWidth(3)
    ps.SetLineColor(col_ps)
    ps.SetMarkerSize(0)
    ps.SetFillStyle(0)

    hs.Add(ps)
    hs.Add(sm)
    #hs.Add(best)

    hs.Draw("nostack hist same")

    bkg.SetFillColor(CreateTransparentColor(12,0.4))
    bkg.SetLineColor(CreateTransparentColor(12,0.4))
    bkg.SetMarkerSize(0)
    bkg.SetMarkerColor(CreateTransparentColor(12,0.4))

    bkg.Draw("e2same")
    data.Draw("E same")


    if mode==1 or True: DrawCMSLogo(pads[0], 'CMS', 'Preliminary', 11, 0.001, -0.07, 0.2, 1.5, '', 1.0)
    else: DrawCMSLogo(pads[0], 'CMS', 'Supplementary', 11, 0.001, -0.07, 0.2, 1.5, '', 1.0)
    DrawTitle(pads[0], '137 fb^{-1} (13 TeV)', 3)

    #Setup legend
    legend = PositionedLegend(0.25,0.3,3,0.02,0.08)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    legend.AddEntry(data,'Data #minus Bkg.',"lep")
    legend.AddEntry(bkg,'Bkg. uncert.',"f")
    legend.AddEntry(sm,'#phi_{#tau#tau} = 0^{#circ}',"l")
    legend.AddEntry(ps,'#phi_{#tau#tau} = 90^{#circ}',"l")
    legend.Draw("same")


    latex.DrawLatex(0.2, 0.85, title)

    c1.SaveAs(plot_name+'.pdf')
