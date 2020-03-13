import os

#year=2018


for year in [2016,2017,2018]:
  out='plot_%(year)i' % vars()
  if year == 2016: extra = '--lumi=\"35.9 fb^{-1} (13 TeV)\" '
  elif year == 2017: extra = '--lumi=\"41.5 fb^{-1} (13 TeV)\" ' 
  else: extra = ''
  for i in range(1,19):
  
    os.system('python scripts/TauIDPlot.py --tauidplot  --file output_shapes/VVLooseVsEle/embed/%(year)i/shapes_%(i)i.root --mode postfit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55 --out=%(out)s_embed --custom_x_range %(extra)s' % vars())
  
    os.system('python scripts/TauIDPlot.py --tauidplot  --file output_shapes/VVLooseVsEle/MC/%(year)i/shapes_%(i)i.root --mode postfit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55 --out=%(out)s_mc --custom_x_range %(extra)s' % vars())
  
    os.system('python scripts/TauIDPlot.py --tauidplot  --file output_shapes/VVLooseVsEle/embed/%(year)i/shapes_%(i)i.root --mode prefit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55 --out=%(out)s_embed --custom_x_range %(extra)s' % vars())
  
    os.system('python scripts/TauIDPlot.py --tauidplot  --file output_shapes/VVLooseVsEle/MC/%(year)i/shapes_%(i)i.root --mode prefit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55 --out=%(out)s_mc --custom_x_range %(extra)s' % vars())
