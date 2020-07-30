from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel
import math

class CPMixture(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.sm_fix = True
        self.do_fa3 = False
        self.useRate = False
        self.kappa = False
        self.ChannelCompatibility = False
        self.VBFangle = False
        self.floatAlpha = False
        self.do2D = False

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("no_sm_fix"):
                self.sm_fix = False
            if po.startswith("do_fa3"):
                self.do_fa3 = True
            if po.startswith("useRate"):
                self.useRate = True
            if po.startswith("ChannelCompatibility"):
                self.ChannelCompatibility = True
            if po.startswith("kappa"):
                self.kappa = True
            if po.startswith("floatAlpha"):
                self.floatAlpha = True
            if po.startswith("do2D"):
                self.do2D = True

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- POI and other parameters ----
          
        params = {
                'pi': math.pi,
                'qqh_sigma1' : 0.2371314,
                'qqh_sigma3' : 2.6707,
                'qqh_sigmamm' : 0.47421,
                'zh_sigma1' : 0.0554268,
                'zh_sigma3' : 2.6708711,
                'zh_sigmamm' : 0.1108536,
                'wh_sigma1' : 0.08607456,
                'wh_sigma3' : 5.6330289,
                'wh_sigmamm' : 0.1721491,
                'a1_qqh_mm' : 1.0,
                'a3_qqh_mm' : 0.297979018705,
                'a1_zh_mm' : 1.0,
                'a3_zh_mm' : 0.144057,
                'a1_wh_mm' : 1.0,
                'a3_wh_mm' : 0.1236136,
                'sm_evt_2016': 1612810.,
                'mm_evt_2016': 1728080.,
                'ps_evt_2016': 1751196.,
                'sm_evt_2017': 11638999.,
                'mm_evt_2017': 1803052.,
                'ps_evt_2017': 1879685.,
                'sigma1_hzz': 1.010599e-2,
                'sigma3_hzz': 1.54779e-3
        }
  
        poiNames = []
 
        if (not self.do_fa3 and not self.ChannelCompatibility and not self.kappa and not self.do2D) or self.floatAlpha:
          self.modelBuilder.doVar('alpha[0,-90,90]') 
          poiNames.append('alpha')
        elif not self.ChannelCompatibility and not self.kappa:
          self.modelBuilder.doVar('fa3[0,0,1]')
          poiNames.append('fa3')

        # if self.kappa:
        #   self.modelBuilder.doVar('kh[1,-2,2]')
        #   poiNames.append('kh')
        #   self.modelBuilder.doVar('ka[0,-2,2]')
        #   poiNames.append('ka')
        #   self.modelBuilder.factory_('expr::muggH("@0*@0 + 2.25*@1*@1", kh, ka)')
        #   if not self.floatAlpha: self.modelBuilder.factory_('expr::alpha("atan(1.5*@1/@0)*180/{pi}", kh, ka)'.format(**params))

        # new kappa scan
        if self.kappa:
            self.modelBuilder.doVar('kappaH[1,-2,2]')
            self.modelBuilder.doVar('kappaA[0,-2,2]')
            poiNames.append('kappaH')
            poiNames.append('kappaA')
            self.modelBuilder.doSet('POI', ','.join(poiNames))
            # scale muggH by kappa^2 + (3/2 * kappa_tilde)^2
            self.modelBuilder.factory_('expr::muggH("@0*@0 + 2.25*@1*@1", kappaH, kappaA)')
            self.modelBuilder.doVar('mutautau[1]')

            self.modelBuilder.factory_('expr::a1("@0", kappaH)'.format(**params))
            self.modelBuilder.factory_('expr::a3("1.5*@0", kappaA)'.format(**params))

            self.modelBuilder.factory_('expr::width_correction("1/(1 + (8.19E-02)*(@0*@0+2.25*@1*@1 -1))", kappaH, kappaA)')

            self.modelBuilder.factory_('expr::sm_scaling("(@0*@0 - @0*@1)*@2", a1, a3, width_correction)')
            self.modelBuilder.factory_('expr::ps_scaling("(@1*@1 - @0*@1)*@2", a1, a3, width_correction)')
            # don't scale MM by 2x SM as XS already scaled by 2
            self.modelBuilder.factory_('expr::mm_scaling("@0*@1*@2", a1, a3, width_correction)') 

        if self.do2D:
          self.modelBuilder.doVar('a2[1,-2,2]')
          poiNames.append('a2')
          self.modelBuilder.doVar('a3[0,-2,2]')
          poiNames.append('a3')
          self.modelBuilder.factory_('expr::muggH("@0*@0 + @1*@1", a2, a3)')

        if self.ChannelCompatibility and not self.kappa and not self.do2D:
          self.modelBuilder.doVar('alpha_em[0,-90,90]')
          poiNames.append('alpha_em')
          self.modelBuilder.doVar('alpha_et[0,-90,90]')
          poiNames.append('alpha_et')
          self.modelBuilder.doVar('alpha_mt[0,-90,90]')
          poiNames.append('alpha_mt')
          self.modelBuilder.doVar('alpha_tt[0,-90,90]')
          poiNames.append('alpha_tt') 

        #self.modelBuilder.doVar('mutautau[1,0,10]')
        self.modelBuilder.doVar('mutautau[1]')
        self.modelBuilder.doVar('muV[1,0,4]')
        if not self.useRate and not self.kappa and not self.do2D:
           self.modelBuilder.doVar('muggH[1,0,10]')
        if not self.VBFangle: self.modelBuilder.doVar('f[0,-1,1]')
        else: self.modelBuilder.doVar('beta[0,-90,90]')

        self.modelBuilder.doSet('POI', ','.join(poiNames))

        self.modelBuilder.doVar('Zero[0]')

        self.modelBuilder.doSet('POI', ','.join(poiNames))


        if self.useRate:
          self.modelBuilder.factory_('expr::muggH("pow(cos(@0/90*{pi}/2),2) + 2.25*pow(sin(@0/90*{pi}/2),2)", alpha)'.format(**params))

        self.modelBuilder.factory_('expr::muV_mutautau("@0*@1", muV, mutautau)')
        self.modelBuilder.factory_('expr::muggH_mutautau("@0*@1", muggH, mutautau)')

        if not self.kappa:
            if not self.ChannelCompatibility:

              if not self.do_fa3:
                if not self.do2D:
                  self.modelBuilder.factory_('expr::a2("sqrt(@0)*cos(@1/90*{pi}/2)", muggH_mutautau, alpha)'.format(**params))
                  self.modelBuilder.factory_('expr::a3("sqrt(@0)*sin(@1/90*{pi}/2)", muggH_mutautau, alpha)'.format(**params))
              else:
                self.modelBuilder.factory_('expr::a3("sqrt( (@0*@1*{sigma1_hzz})/({sigma3_hzz + @1*({sigma1_hzz - sigma3_hzz})}) ) ", muggH_mutautau, fa3)'.format(**params))
                self.modelBuilder.factory_('expr::a2("sqrt(@0-@1*@1)", muggH_mutautau, a3)'.format(**params)) 
  
              self.modelBuilder.factory_('expr::sm_scaling("@0*@0 - @0*@1", a2, a3)')
              self.modelBuilder.factory_('expr::ps_scaling("@1*@1 - @0*@1", a2, a3)')
              self.modelBuilder.factory_('expr::mm_scaling_sync("2*@0*@1", a2, a3)')
              self.modelBuilder.factory_('expr::mm_scaling("@0*@1", a2, a3)')
 
            else: 
              for c in ['em', 'et', 'mt', 'tt']:
                self.modelBuilder.factory_('expr::a2_%s("sqrt(@0)*cos(@1/90*{pi}/2)", muggH_mutautau, alpha_%s)'.format(**params) % (c, c))
                self.modelBuilder.factory_('expr::a3_%s("sqrt(@0)*sin(@1/90*{pi}/2)", muggH_mutautau, alpha_%s)'.format(**params) % (c, c))
   
                self.modelBuilder.factory_('expr::sm_scaling_%s("@0*@0 - @0*@1", a2_%s, a3_%s)' % (c, c, c))
                self.modelBuilder.factory_('expr::ps_scaling_%s("@1*@1 - @0*@1", a2_%s, a3_%s)' % (c, c, c))
                self.modelBuilder.factory_('expr::mm_scaling_%s("@0*@1", a2_%s, a3_%s)' % (c, c, c)) 
                self.modelBuilder.factory_('expr::mm_scaling_sync_%s("2*@0*@1", a2_%s, a3_%s)' % (c, c, c)) 
        
        # For the 0jet and boosted categories since all templates should be the same for SM, MM, and PS sum together all the templates but weight by event numbers in the samples

        self.modelBuilder.doVar('expr::muggH_mutautau_sm_01jetnorm_2016("@0*{sm_evt_2016}/({sm_evt_2016}+{mm_evt_2016}+{ps_evt_2016})", muggH_mutautau)'.format(**params))
        self.modelBuilder.doVar('expr::muggH_mutautau_mm_01jetnorm_2016("@0*{mm_evt_2016}/({sm_evt_2016}+{mm_evt_2016}+{ps_evt_2016})", muggH_mutautau)'.format(**params))
        self.modelBuilder.doVar('expr::muggH_mutautau_ps_01jetnorm_2016("@0*{ps_evt_2016}/({sm_evt_2016}+{mm_evt_2016}+{ps_evt_2016})", muggH_mutautau)'.format(**params))
        self.modelBuilder.doVar('expr::muggH_mutautau_sm_01jetnorm_2017("@0*{sm_evt_2017}/({sm_evt_2017}+{mm_evt_2017}+{ps_evt_2017})", muggH_mutautau)'.format(**params))
        self.modelBuilder.doVar('expr::muggH_mutautau_mm_01jetnorm_2017("@0*{mm_evt_2017}/({sm_evt_2017}+{mm_evt_2017}+{ps_evt_2017})", muggH_mutautau)'.format(**params))
        self.modelBuilder.doVar('expr::muggH_mutautau_ps_01jetnorm_2017("@0*{ps_evt_2017}/({sm_evt_2017}+{mm_evt_2017}+{ps_evt_2017})", muggH_mutautau)'.format(**params))

        # define scaling for JHU VBF/VH templates as a function of muV,mutautau, and fa3 - this uses VBF definition of fa3 

        if not self.VBFangle:
          # dont think this works for the VH not sure why - ? becasue cross sections don't scale in the same way for VH and VBF as is assumed in fa3 paremeterization?
          # i.e sigma_3/sigma_1 not same for WH, VBF and ZH which mean we can have same a1, a3 and muV since muV = a1^2 + a3^2*sigma_3/sigma_1 if sigma_3/sigma_1 changes then soemthing else has to change to compensate
          self.modelBuilder.factory_('expr::signf("@0/abs(@0)", f)')  

          self.modelBuilder.factory_('expr::a1_qqh("sqrt( (1-abs(@0))*@1 )", f, muV_mutautau)')
          self.modelBuilder.factory_('expr::a3_qqh("@0*sqrt( abs(@1)*@2*{qqh_sigma1}/{qqh_sigma3} )", signf, f, muV_mutautau)'.format(**params))

          self.modelBuilder.factory_('expr::vbf_sm_scaling("@0*@0 -@0*@1*{a1_qqh_mm}/{a3_qqh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::vbf_ps_scaling("@1*@1-@0*@1*{a3_qqh_mm}/{a1_qqh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::vbf_mm_scaling("@0*@1/({a1_qqh_mm}*{a3_qqh_mm})", a1_qqh, a3_qqh)'.format(**params)) 

          self.modelBuilder.factory_('expr::a3_zh("@0*sqrt( abs(@1)*@2*{zh_sigma1}/{zh_sigma3} )", signf, f, muV_mutautau)'.format(**params))

          self.modelBuilder.factory_('expr::zh_sm_scaling("@0*@0 -@0*@1*{a1_zh_mm}/{a3_zh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::zh_ps_scaling("@1*@1-@0*@1*{a3_zh_mm}/{a1_zh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::zh_mm_scaling("@0*@1/({a1_zh_mm}*{a3_zh_mm})", a1_qqh, a3_qqh)'.format(**params))

          self.modelBuilder.factory_('expr::wh_sm_scaling("@0*@0 -@0*@1*{a1_wh_mm}/{a3_wh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::wh_ps_scaling("@1*@1-@0*@1*{a3_wh_mm}/{a1_wh_mm}", a1_qqh, a3_qqh)'.format(**params))
          self.modelBuilder.factory_('expr::wh_mm_scaling("@0*@1/({a1_wh_mm}*{a3_wh_mm})", a1_qqh, a3_qqh)'.format(**params))

        else:
          self.modelBuilder.factory_('expr::a1_beta("sqrt(@0)*cos(@1/90*{pi}/2)", muV_mutautau, beta)'.format(**params))
          self.modelBuilder.factory_('expr::a3_beta("sqrt(@0)*sin(@1/90*{pi}/2)", muV_mutautau, beta)'.format(**params))

          self.modelBuilder.factory_('expr::vbf_sm_scaling("@0*@0 - @0*@1", a1_beta, a3_beta)')
          self.modelBuilder.factory_('expr::vbf_ps_scaling("(@1*@1 - @0*@1)*{qqh_sigma1}/{qqh_sigma3}", a1_beta, a3_beta)'.format(**params))
          self.modelBuilder.factory_('expr::vbf_mm_scaling("2*@0*@1*{qqh_sigma1}/{qqh_sigmamm}", a1_beta, a3_beta)'.format(**params))

          self.modelBuilder.factory_('expr::zh_sm_scaling("@0*@0 - @0*@1", a1_beta, a3_beta)')
          self.modelBuilder.factory_('expr::zh_ps_scaling("(@1*@1 - @0*@1)*{zh_sigma1}/{zh_sigma3}", a1_beta, a3_beta)'.format(**params))
          self.modelBuilder.factory_('expr::zh_mm_scaling("2*@0*@1*{zh_sigma1}/{zh_sigmamm}", a1_beta, a3_beta)'.format(**params)) 
       
          self.modelBuilder.factory_('expr::wh_sm_scaling("@0*@0 - @0*@1", a1_beta, a3_beta)')
          self.modelBuilder.factory_('expr::wh_ps_scaling("(@1*@1 - @0*@1)*{wh_sigma1}/{wh_sigma3}", a1_beta, a3_beta)'.format(**params))
          self.modelBuilder.factory_('expr::wh_mm_scaling("2*@0*@1*{wh_sigma1}/{wh_sigmamm}", a1_beta, a3_beta)'.format(**params))


    def getYieldScale(self, bin_, process):

        scalings = []
        if 'ggH' in process and 'hww' not in process:
            chan = ''
            if '_em_' in bin_: chan = 'em'
            elif '_et_' in bin_: chan = 'et'
            elif '_mt_' in bin_: chan = 'mt'
            elif '_tt_' in bin_: chan = 'tt'

            if "ggH_sm" in process or "ggH2jsm" in process or 'reweighted_ggH_htt_0PM' == process:
                scalings.append('sm_scaling')
            elif "ggH_ps" in process or "ggH2jps" in process or 'reweighted_ggH_htt_0M' == process:
                scalings.append('ps_scaling')
            elif "ggH_mm" in process or "ggH2jmm" in process:
                scalings.append('mm_scaling')
            elif 'reweighted_ggH_htt_0Mf05ph0' == process:
                scalings.append('mm_scaling_sync')
            elif 'ggH_ph' in process:
                scalings.append('muggH_mutautau')

            if self.ChannelCompatibility: scalings.append('%s' % chan)

        if ('qqH' in process or 'WH' in process or 'ZH' in process or 'vbf125_powheg' in process or 'wh125_powheg' in process or 'zh125_powheg' in process) and 'hww' not in process:
            if "qqH_sm" in process:
                scalings.append('vbf_sm_scaling')
            elif "qqH_ps" in process:
                scalings.append('vbf_ps_scaling')
            elif "qqH_mm" in process:
                scalings.append('vbf_mm_scaling')
            elif "WH_sm" in process:
                scalings.append('wh_sm_scaling')
            elif "WH_ps" in process:
                scalings.append('wh_ps_scaling')
            elif "WH_mm" in process:
                scalings.append('wh_mm_scaling')
            elif "ZH_sm" in process:
                scalings.append('zh_sm_scaling')
            elif "ZH_ps" in process:
                scalings.append('zh_ps_scaling')
            elif "ZH_mm" in process:
                scalings.append('zh_mm_scaling')
            else:
               scalings.append('muV_mutautau')

        if scalings:
            scaling = '_'.join(scalings)

            if self.sm_fix:
                if "_1_13TeV" in bin_ or "_2_13TeV" in bin_:
                    # simply this by only take SM template for 0 and 1 jet categories
                    if "ggH_sm" in process or 'reweighted_ggH_htt_0PM' in process: scaling = "muggH_mutautau"
                    elif "qqH_htt" in process or 'WH_htt' in process or 'ZH_htt' in process: scaling = "muV_mutautau"
                    else: scaling = "Zero"

                    #year = '2016'
                    #if '2017' in bin_: year = '2017' 
                    #if "ggH_sm" in process: scaling = "muggH_mutautau_sm_01jetnorm_%s" % year
                    #else: scaling = "Zero"
                    #elif "ggH_mm" in process: scaling = "muggH_mutautau_mm_01jetnorm_%s" % year
                    #elif "ggH_ps" in process: scaling = "muggH_mutautau_ps_01jetnorm_%s" % year
                    
                    #if "ggH_sm" in process: scaling = "muggH_mutautau"
                    #elif "ggH_mm" in process: scaling = "Zero" 
                    #elif "ggH_ps" in process: scaling = "Zero"

            if "ggH_ph_htt" in process: scaling = "muggH_mutautau"

            print 'Scaling %s/%s as %s' % (bin_, process,scaling)
            return scaling
        else:
            return 1

CPMixture = CPMixture()


