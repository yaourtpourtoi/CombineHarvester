from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel
import numpy as np

class CPMixtureDecays(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)

    #def setPhysicsOptions(self, physOptions):
    #    for po in physOptions:

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- POI and other parameters ----
            
        poiNames = []
 
        self.modelBuilder.doVar('alpha[0,-1,1]') 
        poiNames.append('alpha')

        self.modelBuilder.doVar('mutautau[1]')
        self.modelBuilder.doVar('muV[1,0,10]')
        self.modelBuilder.doVar('muggH[1,0,10]')

        self.modelBuilder.doSet('POI', ','.join(poiNames))

        params = {
                'pi': np.pi
        }

        self.modelBuilder.factory_('expr::a1("sqrt(@0)*cos(@1/90*{pi}/2)", mutautau, alpha)'.format(**params))
        self.modelBuilder.factory_('expr::a3("sqrt(@0)*sin(@1/90*{pi}/2)", mutautau, alpha)'.format(**params))

        self.modelBuilder.factory_('expr::sm_scaling("@0*@0 - @0*@1", a1, a3)')
        self.modelBuilder.factory_('expr::ps_scaling("@1*@1 - @0*@1", a1, a3)')
        self.modelBuilder.factory_('expr::mm_scaling("2*@0*@1", a1, a3)')

        self.modelBuilder.factory_('expr::muV_mutautau("@0*@1", muV, mutautau)')
        for x in ['sm_scaling', 'ps_scaling', 'mm_scaling']:
          self.modelBuilder.factory_('expr::muV_%s("@0*@1", muV, %s)' % (x,x))
          self.modelBuilder.factory_('expr::muggH_%s("@0*@1", muggH, %s)' % (x,x))

    def getYieldScale(self, bin_, process):

        scalings = []

        if 'ggH' in process: 
            scalings.append('muggH')
        if 'qqH' in process or 'ZH' in process or 'WH' in process: 
            scalings.append('muV')
        
        if 'ZH' in process or 'WH' in process:
            scalings.append('mutautau')

        if ('ggH' in process or 'qqH' in process) and 'hww' not in process:
            if "sm" in process:
                scalings.append('sm_scaling')
            elif "ps" in process:
                scalings.append('ps_scaling')
            elif "mm" in process:
                scalings.append('mm_scaling')

        if scalings:
            scaling = '_'.join(scalings)

            print 'Scaling %s/%s as %s' % (bin_, process,scaling)
            return scaling
        else:
            return 1

CPMixtureDecays = CPMixtureDecays()


