import os

for i in range(1,11):

  os.system('python scripts/TauIDPlot.py --tauidplot  --file tauid_shapes/shapes_embed_%(i)i.root --mode postfit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55' % vars())

  os.system('python scripts/TauIDPlot.py --tauidplot  --file tauid_shapes/shapes_MC_%(i)i.root --mode postfit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55' % vars())

  os.system('python scripts/TauIDPlot.py --tauidplot  --file tauid_shapes/shapes_embed_%(i)i.root --mode prefit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55' % vars())

  os.system('python scripts/TauIDPlot.py --tauidplot  --file tauid_shapes/shapes_MC_%(i)i.root --mode prefit --x_title \"m_{vis} (GeV)\" --ratio --channel=mt --extra_pad=0.55' % vars())
