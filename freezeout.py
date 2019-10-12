# freeze-out implementation in python
# latest update 20190923 1843
from __future__ import print_function

import numpy as np
from numpy import *
import scipy.integrate
from math import sqrt, log, atan2, pi
from cmath import exp
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import os, sys
import traceback
import glob

import threading
import gi
gi.require_version('Gtk', '2.0')
from gi.repository import Gtk, Gdk, GdkPixbuf

iverbose = 4

class nameConventions:

    def __init__(self,  name='auau-nexus.fzo' ):
        self.name = name
        self.hbarc = 1.97
        self.c = 1.0

    def makeName( self, i = 1 ):
        return self.name.replace('.fzo', '-%d.fzo'%i )

    def outputFileName( self, label1, label2='', label3='' ):
        labelme = str(label1)
        if label2 != '':
            labelme += '_'+str(label2)
            if label3 != '':
                labelme += '_'+str(label3)
        return 'zzout_'+str(labelme)+'.dat'

class FzoDramaFactory:

    def __init__( self ):
        self.eventCounter = 0
        self.fileIsMissing = True
        name = 'x'
        if self.islocked():
            print('WARNING: file %s is locked, we will remove the lock file(s) and proceed'%name)
            self.cleanLocks( os.getcwd() )
        self.lock()

    def __del__( self ):
        FzoDramaFactory.cleanLocks( os.getcwd() )

    def lock(self):
        open( 'x'+'.lock', 'w')

    def islocked(self):
        return os.path.exists( 'x'+'.lock' )

    @staticmethod
    def cleanLocks( dir ):
        print( glob.glob( dir+'/*.lock' ) )
        for lock in glob.glob( dir+'/*.lock' ):
            print('removing lock', lock)
            os.remove( lock )

class FreeZeOutSPH: #a frozen-out SPH particle which store the information of a small element of freeze-out surface
    def __init__( self, randomized=False, in_tauf=None,in_xf=None,in_yf=None,in_etaf=None,in_d3SigmaTau=None,in_d3SigmaX=None,in_d3SigmaY=None,in_d3SigmaEta=None,in_uTau=None,in_uX=None,in_uY=None,in_uEta=None,in_ef=None,in_pf=None,in_Tf=None,in_muBf=None,in_muSf=None): #tauf , xf , yf , etaf , d3SigmaTau , d3SigmaX, d3SigmaY, d3SigmaEta, uTau , uX, uY, uEta, ef , pf, Tf , muBf, components of viscosities
        saved_args = locals()
        if randomized == True:
            self.tauf = 1.0
            self.xf = 0.0
            self.yf=0.0
            self.etaf=0.0
            self.d3SigmaTau=5.0+random.uniform(0,1)
            self.d3SigmaX=0.0+random.uniform(-1,1)
            self.d3SigmaY=0.0+random.uniform(-1,1)
            self.d3SigmaEta=0.0+random.uniform(-1,1)
            self.uX=0.0+random.uniform(-1,1)*0.2
            self.uY=0.0+random.uniform(-1,1)*0.2
            self.uTau=sqrt(1.0+self.uX**2+self.uY**2)+random.uniform(0,1)*0.2
            self.uEta=sqrt(self.uTau**2-self.uX**2-self.uY**2-1.0)/self.tauf
            self.ef=1.0+random.uniform(0,1)*0.5
            self.pf=self.ef/3.0
            self.Tf=0.15
            self.muBf=0.0
            self.muSf=0.0
        elif str(saved_args).count('None') != 0 :
            print('fatal error in FreeZeOutSPH.__init__(), program halt!')
            exit(1)
        else:
            self.tauf = in_tauf
            self.xf = in_xf
            self.yf= in_yf
            self.etaf= in_etaf
            self.d3SigmaTau=in_d3SigmaTau
            self.d3SigmaX=in_d3SigmaX
            self.d3SigmaY=in_d3SigmaY
            self.d3SigmaEta=in_d3SigmaEta
            self.uTau=in_uTau
            self.uX=in_uX
            self.uY=in_uY
            self.uEta=in_uEta
            self.ef=in_ef
            self.pf=in_pf
            self.Tf=in_Tf
            self.muBf=in_muBf
            if in_muSf is not None:
                self.muSf=in_muSf
            else:
                self.muSf=0.0

class FreeZeOutSurface: #vessle class to store freeze-out surface elements
    def __init__( self, opt='Fake' ):
        print ('FreezeOutSurface Class instantiated')
        self.FZOList = []
        self.numFZOSPH = 0
        if opt=='Fake':
            print('randomly creating fake freeze-out surface...')
            for ifzsph in range(50):
                fzsph = FreeZeOutSPH( randomized=True )
                self.addFZOSPH(fzsph)
            print('done! ',self.numFZOSPH,' freeze-out sph particles have been ceated!')

    def addFZOSPH( self, anFZOSPH=None ):
        if anFZOSPH is not None:
            if isinstance(anFZOSPH, FreeZeOutSPH):
                self.FZOList.append(anFZOSPH)
                self.numFZOSPH += 1
            else:
                print ('fatal type error in FreeZeOutSurface.addFZOSPH(), program halt!')
                exit(1)
        else:
            print ('critical warning in FreeZeOutSurface.addFZOSPH(), no input found!')
    def FZOSPH(self, indexSPH=None ):
        if indexSPH is not None:
            return self.FZOList[indexSPH]

class ElementaryParticle:
    def __init__( self, m=None, s=None, q=None, b=None, g=None, sgl=None ):
        print ('ElementaryParticle Class instantiated')
        saved_args = locals()
        if str(saved_args).count('None') == 0 :
            self.m = m #mass
            self.s = s #strangeness
            self.q = q #electric charge
            self.b = b #baryon number
            self.g = g #degeneracy
            self.sgl = sgl #boson (+1) or fermion (-1)
        else:
            print('fatal error in ElementaryParticle.__init__(), program halt!')
            exit(1)

class ParticleDataTable:
    def __init__( self, opt='Default' ):
        print ('ParticleDataTable Class instantiated')
        self.pIDList = []
        self.numIDList = 0
        if opt=='Default':
            pion0=ElementaryParticle(m=0.135, s=0.0, q=0.0, b=0.0, g=1.0, sgl=1)
            self.addElementaryParticle(anElementaryParticle=pion0)
    def addElementaryParticle( self, anElementaryParticle=None ):
        if anElementaryParticle is not None:
            if isinstance(anElementaryParticle, ElementaryParticle):
                self.pIDList.append(anElementaryParticle)
                self.numIDList += 1
            else:
                print('fatal type error in ParticleDataTable.anElementaryParticle(), program halt!')
                exit(1)
        else:
            print('critical warning from ParticleDataTable.anElementaryParticle(), no input found!')
        return self
    def EleParticle( self, ipindex=None ):
        return self.pIDList[ipindex]

class ParticleOutputList:
    def __init__( self, opt='Default' ):
        print ('ParticleOutputList Class instantiated')
        self.pOutputID = []
        self.numPOutputList = 0
        if opt=='Default':
            self.addOutputParticle(0) #add pion
    def addOutputParticle( self, anOutputParticleID=None ):
        if anOutputParticleID is not None:
            if isinstance(anOutputParticleID, int):
                self.pOutputID.append(anOutputParticleID)
                self.numPOutputList += 1
            else:
                print('fatal type error in ParticleOutputList.anOutputParticleID(), program halt!')
                exit(1)
        else:
            print('critial warning from ParticleOutputList.addOutputParticle(), no input found!')
        return self
    def iOutParticleID( self, iOutpindex=None ):
        return self.pOutputID[iOutpindex]

class SpectrumFlowEvaluation:

    def __init__( self, FZOsrfc='Fake', PDTable='Default', POutList='Default', nEtaBin=None, etaMax=None, etaMin=None, nYBin=None, yMax=None, yMin=None, nPtBin=None, ptMax=None, ptMin=None, nPhiBin=None, nHarmonics=None): #python 2.7
        print ('SpectrumFlowEvaluation Class instantiated')
        self.doneHeader = False
        self.gridsHaveBeenDefined = False
        if FZOsrfc =='Fake':
            FZOsrfc = FreeZeOutSurface('Fake')
        if PDTable == 'Default':
            PDTable = ParticleDataTable('Default')
        if POutList == 'Default':
            POutList = ParticleOutputList('Default')
        if not (isinstance(FZOsrfc, FreeZeOutSurface) and isinstance(PDTable, ParticleDataTable) and isinstance(POutList, ParticleOutputList)):
            print ('Fatal error in SpectrumFreezeOut initialization, program halt!')
            exit(1)
        if self.gridsHaveBeenDefined == False: #first time runner
            saved_args = locals()
            self.FZOsrfc = FZOsrfc
            self.PDTable = PDTable
            self.PoutList = POutList
            self.nEtaBin = 120
            self.etaMax = 7.0
            self.etaMin = -self.etaMax
            self.nYBin = 120
            self.yMax = 7.0
            self.yMin = -self.yMax
            self.ptMax = 5.0
            self.ptMin = 0.0
            self.nPtBin = 12
            self.phiMax = 2.0*np.pi
            self.phiMin = 0.0
            self.nPhiBin = 12
            self.nHarmonics = 4
            if str(saved_args).count('None') != 11 :print('Spectrum grids dimension has been modified:')
            if nEtaBin is not None:
                self.nEtaBin = int(nEtaBin)
                print('nEtaBin:',self.nEtaBin)
            if etaMax is not None:
                self.etaMax = float(etaMax)
                print('etaMax:',self.etaMax)
            if etaMin is not None:
                self.etaMin = float(etaMin)
                print('etaMin:',self.etaMin)
            if nYBin is not None:
                self.nYBin = int(nYBin)
                print('nEtaBin:',self.nEtaBin)
            if yMax is not None:
                self.yMax = float(yMax)
                print('etaMax:',self.etaMax)
            if yMin is not None:
                self.yMin = float(yMin)
                print('etaMin:',self.etaMin)
            if nPtBin is not None:
                self.nPtBin = int(nPtBin)
                print('nPtBin:',self.nPtBin)
            if ptMax is not None:
                self.ptMax = float(ptMax)
                print('ptMax:',self.ptMax)
            if ptMin is not None:
                self.ptMin = float(ptMin)
                print('ptMin:',self.ptMin)
            if nPhiBin is not None:
                self.nPhiBin = int(nPhiBin)
                print('nPhiBin:',self.nPhiBin)
            if nHarmonics is not None:
                self.nHarmonics = int(nHarmonics)
            self.deltaEta = (self.etaMax-self.etaMin)/self.nEtaBin
            self.deltaY = (self.yMax-self.yMin)/self.nYBin
            self.deltaPt = (self.ptMax-self.ptMin)/self.nPtBin
            self.deltaPhi = (self.phiMax-self.phiMin)/self.nPhiBin
            self.D3nDyptDptDphi = zeros([self.PoutList.numPOutputList,self.nPtBin, self.nYBin, self.nPhiBin],dtype=np.float32) #np.float64 is double precision
            self.D3nDetaptDptDphi = zeros([self.PoutList.numPOutputList,self.nPtBin, self.nYBin, self.nPhiBin],dtype=np.float32)
            self.VnYPtEP = zeros([self.nHarmonics, self.PoutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            self.VnYEtaEP = zeros([self.nHarmonics, self.PoutList.numPOutputList+1, self.nEtaBin],dtype=np.float32)
            self.VnYPtCMLT = zeros([self.nHarmonics, self.PoutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            self.VnYEtaCMLT = zeros([self.nHarmonics, self.PoutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            self.gridsHaveBeenDefined = True
            print('grids definition:')
            print('nEtaBin: ',self.nEtaBin,' etaMax: ',self.etaMax, ' etaMin: ',self.etaMin)
            print('nYBin: ',self.nYBin,' yMax: ',self.yMax,' yMin: ',self.yMin)
            print('nPtBin: ',self.nPtBin,' ptMax: ',self.ptMax,' ptMin: ',self.ptMin)
            print('nPhiBin: ',self.nPhiBin,' phiMax: ',self.phiMax,' phiMin: ',self.phiMin)
            print('nHarmonics: ',self.nHarmonics)

    def header( self, nameset ):
        f = open( nameset.outputFileName('pt'), 'w' )
        f.write('pt'+' '+'vn'+'\n')
        f.close()
        f = open( nameset.outputFileName('eta'), 'w' )
        f.write('eta'+' '+'vn'+'\n')
        f.close()
        self.doneHeader = True

    def SpectrumFlowWriteOutputFile( self, nameset ):
        print ('SpectrumFlowOutput Class instantiated')
        if not self.doneHeader:
            self.header( nameset )
        f = open( nameset.outputFileName('pt'), 'a' ) #'a' for append
        for iPtBin in range(self.nPtBin):
            pt = self.ptMin+(iPtBin+0.5)*(self.ptMax-self.ptMin)
            f.write(str(pt)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYPtEP[:, self.PoutList.numPOutputList, iPtBin]])+'\n')
            if iverbose > 8: print(str(pt)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYPtEP[:, self.PoutList.numPOutputList, iPtBin]])+'\n')
        f.close()
        f = open( nameset.outputFileName('eta'), 'a' ) #'a' for append
        for iEtaBin in range(self.nEtaBin):
            eta = self.etaMin+(iEtaBin+0.5)*(self.etaMax-self.etaMin)
            f.write(str(eta)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYEtaEP[:, self.PoutList.numPOutputList, iEtaBin]])+'\n')
            if iverbose > 8: print(str(eta)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYEtaEP[:, self.PoutList.numPOutputList, iEtaBin]])+'\n')
        f.close()
        print('The obtained results have been saved to the file(s).')
        return self

    def evaluateGridD3nDyptDptDphi( self, nameset ):
        print ('Subroutine evaluateGridD3nDyptDptDphi() invoked')
        for iescolTmp in range(self.FZOsrfc.numFZOSPH):
            for iParticleSpeciesOutput in range(self.PoutList.numPOutputList):
                for iPtBin in range(self.nPtBin):
                    rptTmp=self.ptMin+(iPtBin+0.5)*(self.ptMax-self.ptMin)/float(self.nPtBin)
                    for iYBin in range(self.nYBin):
                        ryTmp=self.yMin+(iYBin+0.5)*(self.yMax-self.yMin)/float(self.nYBin)
                        for iPhiBin in range(self.nPhiBin):
                            rphiTmp=self.phiMin+(iPhiBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiBin)
                            tmpDis=self.fD3nDptDyDphi(nameset,rpt=rptTmp,ry=ryTmp,rphi=rphiTmp,iescol=iescolTmp,ipid=iParticleSpeciesOutput)
                            self.D3nDyptDptDphi[iParticleSpeciesOutput,iPtBin,iYBin,iPhiBin]+=tmpDis
        return self

    def evaluateOutputD3nDyptDptVnVsPt( self ):
        print ('Subroutine evaluateOutputD3nDyptDptVnVsPt() invoked')
#        psib_temp     = psib
#        psia_temp     = psia
#        phievent_temp = phievent
        psib          = 0.0
        psia          = 0.0
        phievent      = 0.0
        sum0          = 0.0
        sum0pt        = 0.0
        phievt        = zeros([self.nHarmonics],dtype=np.float32)
        psisin        = zeros([self.nHarmonics],dtype=np.float32)
        psicos        = zeros([self.nHarmonics],dtype=np.float32)
        psisinwpt     = zeros([self.nHarmonics],dtype=np.float32)
        psicoswpt     = zeros([self.nHarmonics],dtype=np.float32)
        for iPtBin in range(self.nPtBin):
            sum0a5   = 0.0
            sum0a5v1 = 0.0
            sum0a5v2 = 0.0
            soma0a56 = 0.0
            soma0a57 = 0.0
            for iHarmonics in range(self.nHarmonics):
                soma0a5cos = zeros([self.nHarmonics],dtype=np.float32)
            for iParticleSpeciesOutput in range(self.PoutList.numPOutputList):
                rpt=self.ptMin+(iPtBin+0.5)*(self.ptMax-self.ptMin)/float(self.nPtBin)
                if (iPtBin == 0 or iPtBin == self.nPtBin-1):
                    fator0 = 1.0
                elif iPtBin%2 == 0:
                    fator0 = 4.0
                elif iPtBin%2 == 1:
                    fator0 = 2.0
                sum1  = 0.0
                sum1pt= 0.0
                sum1v1= 0.0
                sum1v2= 0.0
                soma16= 0.0
                soma17= 0.0
                for iHarmonics in range(self.nHarmonics):
                    soma1sin    = zeros([self.nHarmonics],dtype=np.float32)
                    soma1cos    = zeros([self.nHarmonics],dtype=np.float32)
                    soma1sinwpt = zeros([self.nHarmonics],dtype=np.float32)
                    soma1coswpt = zeros([self.nHarmonics],dtype=np.float32)
                for iPhiBin in range(self.nPhiBin):
                    rphi=self.phiMin+(iPhiBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiBin)
                    fator1=3.0 #periodic function
                    sum2 = 0.0
                    sum2pt=0.0
                    for iHarmonics in range(self.nHarmonics):
                        soma2sin    = zeros([self.nHarmonics],dtype=np.float32)
                        soma2cos    = zeros([self.nHarmonics],dtype=np.float32)
                        soma2sinwpt = zeros([self.nHarmonics],dtype=np.float32)
                        soma2coswpt = zeros([self.nHarmonics],dtype=np.float32)
                        for iYBin in range(self.nYBin):
                            ry=self.yMin+(iYBin+0.5)*(self.yMax-self.yMin)/float(self.nYBin)
                            if (ry < 0.0):
                                fase=psib
                            else:
                                fase=psia
                            if (iYBin == 1 or iYBin == self.nYBin-1):
                                fator2 = 1.0
                            elif iYBin%2 == 0:
                                fator2 = 4.0
                            elif iYBin%2 == 1:
                                fator2 = 2.0
                            tmp2    = self.D3nDyptDptDphi[iParticleSpeciesOutput,iPtBin,iYBin,iPhiBin]
                            sum2   += tmp2*fator2
                            sum2pt += tmp2*rpt*fator2
                            for iHarmonics in range(self.nHarmonics):
                                soma2sin[iHarmonics]+=tmp2*fator2*sin(float(iHarmonics)*(rphi-phievt[iHarmonics]))
                                soma2cos[iHarmonics]+=tmp2*fator2*cos(float(iHarmonics)*(rphi-phievt[iHarmonics]))
                                soma2sinwpt[iHarmonics]+=tmp2*rpt*fator2*sin(float(iHarmonics)*(rphi-phievt[iHarmonics]))
                                soma2coswpt[iHarmonics]+=tmp2*rpt*fator2*cos(float(iHarmonics)*(rphi-phievt[iHarmonics]))

                        sum1   += sum2*fator1
                        sum1pt += sum2pt*fator1
                        sum1v1 += sum2*fator1*cos(rphi-fase)
                        sum1v2 += sum2*fator1*cos(2.0*(rphi-fase))
                        soma16 += sum2*fator1*cos(2.0*(rphi-phievent))
                        soma17 += sum2*fator1*cos(rphi-phievent)
                        for iHarmonics in range(self.nHarmonics):
                            soma1sin[iHarmonics]    += soma2sin[iHarmonics]*fator1
                            soma1cos[iHarmonics]    += soma2cos[iHarmonics]*fator1
                            soma1sinwpt[iHarmonics] += soma2sinwpt[iHarmonics]*fator1
                            soma1coswpt[iHarmonics] += soma2coswpt[iHarmonics]*fator1
                    sum1   *= self.deltaY/3.0
                    sum1   *= self.deltaPhi/3.0
                    sum1pt *= self.deltaY/3.0
                    sum1pt *= self.deltaPhi/3.0
                    sum1v1 *= self.deltaY/3.0
                    sum1v1 *= self.deltaPhi/3.0
                    sum1v2 *= self.deltaY/3.0
                    sum1v2 *= self.deltaPhi/3.0
                    soma16 *= self.deltaY/3.0
                    soma16 *= self.deltaPhi/3.0
                    soma17 *= self.deltaY/3.0
                    soma17 *= self.deltaPhi/3.0
                    for iHarmonics in range(self.nHarmonics):
                        soma1sin[iHarmonics] *= self.deltaY/3.0
                        soma1sin[iHarmonics] *= self.deltaPhi/3.0
                        soma1cos[iHarmonics] *= self.deltaY/3.0
                        soma1cos[iHarmonics] *= self.deltaPhi/3.0
                sum0a5    += sum1 #only to sum up different specials, therefore no fator0
                sum0a5v1  += sum1v1
                sum0a5v2  += sum1v2
                soma0a56  += soma16
                soma0a57  += soma17

                sum0      += sum1*fator0 #this summation should take place for every species
                sum0pt    += sum1pt*fator0
                for iHarmonics in range(self.nHarmonics):
                    soma0a5cos[iHarmonics] += soma1cos[iHarmonics]
                    psisin[iHarmonics]     += soma1sin[iHarmonics]*fator0
                    psicos[iHarmonics]     += soma1cos[iHarmonics]*fator0
                    psisinwpt[iHarmonics]  += soma1sinwpt[iHarmonics]*fator0
                    psicoswpt[iHarmonics]  += soma1coswpt[iHarmonics]*fator0

                self.VnYPtEP[0,iParticleSpeciesOutput,iPtBin]       += sum1 #multiplicity
                self.VnYPtEP[0,self.PoutList.numPOutputList,iPtBin] += sum1 #all particle species
                self.VnYPtEP[1,iParticleSpeciesOutput,iPtBin]       += sum1v1/sum1 #v1
                self.VnYPtEP[2,iParticleSpeciesOutput,iPtBin]       += sum1v2/sum1 #v2

            self.VnYPtEP[1,self.PoutList.numPOutputList,iPtBin] += sum0a5v1/sum0a5 #v1, all particle species
            self.VnYPtEP[2,self.PoutList.numPOutputList,iPtBin] += sum0a5v2/sum0a5 #v2, all particle species

        sum0   *= self.deltaPt/3.0
        sum0pt *= self.deltaPt/3.0
        for iHarmonics in range(self.nHarmonics):
            psisin[iHarmonics] *= self.deltaPt/3.0
            psicos[iHarmonics] *= self.deltaPt/3.0
        return self

    def fD3nDptDyDphi(self,nameset,rpt=None,ry=None,rphi=None,iescol=None,ipid=None):
        saved_args = locals()
        if str(saved_args).count('None') != 0 :
            print('fD3nDyDpTDeta() input fatal error, program halt!')
            exit(1)

        detetmin = 0.170             #the highest possible value
        Tfzo = min(detetmin,self.FZOsrfc.FZOSPH(iescol).Tf)
        ub = self.FZOsrfc.FZOSPH(iescol).muBf
        us = self.FZOsrfc.FZOSPH(iescol).muSf
        b  = self.PDTable.EleParticle(ipid).b
        m  = self.PDTable.EleParticle(ipid).m
        g  = self.PDTable.EleParticle(ipid).g
        ssq  = self.PDTable.EleParticle(ipid).s
        sinal = self.PDTable.EleParticle(ipid).sgl
        uu = b*ub+ssq*us
        mt = sqrt(rpt**2+m**2)

# 4-vetor s_{\mu} ({\tau}xy{\eta}-coord)
        stau   = self.FZOsrfc.FZOSPH(iescol).d3SigmaTau
        sx     = self.FZOsrfc.FZOSPH(iescol).d3SigmaX
        sy     = self.FZOsrfc.FZOSPH(iescol).d3SigmaY
        seta   = self.FZOsrfc.FZOSPH(iescol).d3SigmaEta
# 4-vetor u_{\mu} ({\tau}xy{\eta}-coord)
        utau   = self.FZOsrfc.FZOSPH(iescol).uTau
        ux     = self.FZOsrfc.FZOSPH(iescol).uX
        uy     = self.FZOsrfc.FZOSPH(iescol).uY
        ueta   = self.FZOsrfc.FZOSPH(iescol).uEta

        tau    = self.FZOsrfc.FZOSPH(iescol).tauf
        x      = self.FZOsrfc.FZOSPH(iescol).xf
        y      = self.FZOsrfc.FZOSPH(iescol).yf
        eta    = self.FZOsrfc.FZOSPH(iescol).etaf
# 4-vetor st_{\mu} (txyz-coord)
        st     =  cosh(eta)*stau-1.0/tau*sinh(eta)*seta
        sx     =  sx
        sy     =  sy
        sz     = -sinh(eta)*stau+1.0/tau*cosh(eta)*seta
# 4-vetor u^{\mu} (txyz-coord)
        ut     =  cosh(eta)*utau-1.0/tau*sinh(eta)*ueta
        ux     = -ux
        uy     = -uy
        uz     =  sinh(eta)*utau-1.0/tau*cosh(eta)*ueta

        vx     = ux/ut
        vy     = uy/ut
        vz     = uz/ut

        v2     = vx*vx+vy*vy+vz*vz
        gamma  = 1.0/sqrt(1.0-v2)
        if (1.0-v2) < 0.0:
            print('NaN in d3n! program halt!')
            print('(1.d0-v2)=',(1.0-v2))
            print('utau,ux,uy,ueta=',utau,ux,uy,ueta)
            print('tau,x,y,eta=',tau,x,y,eta)
            exit(1)

# 4-vetor p^{\mu} (txyz-coord); lab frame
#       p0=sqrt((rpt*cosh(reta))**2+m**2)
        p0=mt*cosh(ry)
        px=rpt*cos(rphi)
        py=rpt*sin(rphi)
#       pz=rpt*sinh(reta)
        pz=mt*sinh(ry)

        pu=gamma*(p0-px*vx-py*vy-pz*vz)

        fex   = exp((-pu+uu)/Tfzo)
        distr = fex/(1.0+sinal*fex)
        if (((-pu+uu)/Tfzo)>(888.0)): #python sys.float_info, max_exp=1024
            distr=1.0
        elif (((-pu+uu)/Tfzo)<(-888.0)):
            distr=0.0

        taa=0.5 * (sign(st*p0+sx*px+sy*py+sz*pz) + 1.0)
        taa=taa * (st*p0+sx*px+sy*py+sz*pz)
#        taa=taa*(rpt**2*cosh(reta))
#        taa=taa/sqrt(rpt**2*cosh(reta)**2+m**2)    #see SPheRIO remarks and key notes on hydro
        taa=taa*rpt
        taa=taa*g/(2.0*pi)**3*distr
        taa=taa/(nameset.hbarc**3)
        if taa.imag != 0.0:
            print('critial warning from fD3nDptDyDphi(), result is not real!')
        d3ndptdydphi=taa.real
        return d3ndptdydphi

    def evaluateOutputD3nDetaDpTDeta():
        print ('Subroutine evaluateOutputDn3DetapTdpT() invoked')
        pass
        return self

    def evaluateOutputDnDeta():
        print ('Subroutine evaluateOutputDnDeta() invoked')
        pass
        return self

class ParticleDecay:
    def __init__( self ):
        print ('ParticleDecay Class instantiated')

class HadronizationMonteCarlo:
    def __init__( self ):
        print ('HadronizationMonteCarlo Class instantiated')

class PlotOnScreen:
    def __init__( self ):
        print ('PlotOnScreen Class instantiated')

    def readOutputFile( self, nameset ):
        pt = []
        dnptdpt = []
        v1pt = []
        v2pt = []
        eta = []
        dndeta = []
        v1eta = []
        v2eta = []

        data = open( nameset.outputFileName('pt'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ') # reorganize the obtained line into strings separated by single blank ' '
            if linesplit[0] == 'pt':
                lambdaPt= 1
                continue
            else:
                pt += [ float(linesplit[0]) ]
                dnptdpt += [ float(linesplit[1]) ]
                v1pt += [ float(linesplit[2]) ]
                v2pt += [ float(linesplit[3]) ]
        data = open( nameset.outputFileName('eta'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
            if linesplit[0] == 'eta':
                lambdaEta= 2
                continue
            else:
                eta += [ float(linesplit[0]) ]
                dndeta += [ float(linesplit[1]) ]
                v1eta += [ float(linesplit[2]) ]
                v2eta += [ float(linesplit[3]) ]
        pt = array(pt)
        dnptdpt = array(dnptdpt)
        v1pt = array(v1pt)
        v2pt = array(v2pt)
        eta = array(eta)
        dndeta = array(dndeta)
        v1eta = array(v1eta)
        v2eta = array(v2eta)
        if iverbose > 8: print('pt, dnptdpt, v1pt, v2pt\n',pt, dnptdpt, v1pt, v2pt)
        return pt, dnptdpt, pt, v2pt, eta, dndeta, eta, v2eta

    def plotOutputFile( self, nameset  ):
        for plotType in ["pt-spectrum","eta-spectrum"]:
            plt.rcParams['font.size'] = 12.
            plt.rcParams['font.family'] = "serif"
            plt.rcParams["xtick.labelsize"] = 'xx-small'
            plt.rcParams["ytick.labelsize"] = 'xx-small'
            plt.rcParams["ytick.labelsize"] = 'xx-small'
            plt.rcParams['text.usetex'] = True
            plt.rcParams['text.latex.unicode'] = True
            width = plt.rcParams['figure.figsize'][0]/2

            fig = plt.figure( figsize = (3*width,width*.75) )
            grid = plt.GridSpec(2, 2, hspace=0, wspace=0)
            subFig1 = fig.add_subplot(grid[0,0])
            subFig2 = fig.add_subplot(grid[1,0], sharex=subFig1)
            subFig3 = fig.add_subplot(grid[0,1])
            subFig4 = fig.add_subplot(grid[1,1], sharex=subFig3)
            fig.subplots_adjust(wspace=0)
            subFig1.grid(True)
            subFig2.grid(True)
            subFig3.grid(True)
            subFig4.grid(True)

            x1, y1, x2, y2, x3, y3, x4, y4 = self.readOutputFile( nameset )
            labelMe = 'v'
            indexMe = 2

            subFig1.plot(x1[:], y1[:], '-b', lw=.5, label=r'$n$={:f}'.format(indexMe) )
            subFig2.plot(x2[:], y2[:], '-b', lw=.5, label=' ' )
            subFig3.plot(x3[:], y3[:], '-r', lw=.5,  label=r'$n$={:f}'.format(indexMe) )
            subFig4.plot(x4[:], y4[:], '-r', lw=.5,  label=' ' )

            subFig1.set_ylabel(r'$v_n$', fontsize=16)
            subFig2.set_ylabel(r'$v_n$', fontsize=16)
            subFig2.set_xlabel(r'$p_T$', fontsize=16)
            subFig4.set_xlabel(r'$\eta$', fontsize=16)

            subFig1.locator_params(nbins=4, axis='y')
            subFig2.locator_params(nbins=4, axis='y')
            subFig3.locator_params(nbins=4, axis='y')
            subFig4.locator_params(nbins=4, axis='y')

            subFig1.legend(loc='lower right')
            subFig2.legend(loc='lower right')
            subFig3.legend(loc='lower right')
            subFig4.legend(loc='lower right')

            subFig1.legend( fancybox=True, framealpha=0 )
            subFig2.legend( fancybox=True, framealpha=0 )
            subFig3.legend( fancybox=True, framealpha=0 )
            subFig4.legend( fancybox=True, framealpha=0 )

            plt.setp(subFig2.get_yticklabels(), visible=True)
            plt.setp(subFig1.get_xticklabels(), visible=False)
            plt.setp(subFig4.get_yticklabels(), visible=True)
            plt.setp(subFig3.get_xticklabels(), visible=False)

            fig.tight_layout()
            if not os.path.exists('figures/'): os.makedirs('figures/')
            outname = 'figures/f_%s%s.dat'%(labelMe,indexMe)
            figname = outname.replace('.dat', '.png')
            fig.savefig( figname )
            print('plot saved to', figname )
            plt.close()
        print('The obtained results have been plotted to the file(s).')
        return self

    def showPlots(self):
        if iverbose > 2: print ('*showFlowPlot() invoked!')
        app = drawFlowPlotApp()
        Gtk.main()

class drawFlowPlotApp(Gtk.Window):

    def __init__(self):
        Gtk.Window.__init__( self, title="PCA plots" )
        self.connect( 'destroy', Gtk.main_quit )
        self.set_border_width(3)
        self.maximize()

        self.vbox = Gtk.VBox()

        self.parameterPanels = []
        self.paramsList = ['type', 'index' ] #type of plot, harmonic index
        self.paramsEntry = []
        self.images = []
        self.threads = []
        self.newParameterPanel( self.parameterPanels )
        self.newParameterPanel( self.parameterPanels )
        self.newParameterPanel( self.parameterPanels, empty = True )

        scrolled = Gtk.ScrolledWindow()
        scrolled.add_with_viewport( self.vbox )
        scrolled.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        self.add( scrolled )

        self.show_all()

    def onaddclicked( self, widget ):
        self.parameterPanels[-1].destroy()
        del self.parameterPanels[-1]
        self.newParameterPanel( self.parameterPanels )
        self.newParameterPanel( self.parameterPanels, empty = True )
        self.show_all()

    def onminusclicked( self, widget, box ):
        box.destroy()
        ind = self.parameterPanels.index(box)
        self.parameterPanels[ ind ].destroy()
        del self.parameterPanels[ ind ]
        del self.paramsEntry[ ind ]
        self.show_all()

    def newParameterPanel(self, panel, empty = False):
        box = Gtk.HBox()
        panel.append( box )
        this = panel[-1]
        self.vbox.pack_start(box, False, False, 3 )

        if empty:
            button = Gtk.Button()
            box.pack_start( button,  False, False, 3 )
            button.set_label(' + ')
            button.connect( 'clicked', self.onaddclicked )
            return

        options = Gtk.VBox()
        box.pack_start( options, False, False, 3 )

        firstLine = Gtk.HBox()
        title = Gtk.Label()
        destroyButton = Gtk.Button()
        options.pack_start(firstLine, False, False, 3 )

        firstLine.pack_start(title, True, True, 3)
        firstLine.pack_end(destroyButton, False, False, 3)

        title.set_label('parameters')
        destroyButton.set_label(' x ')
        destroyButton.connect( 'clicked', self.onminusclicked, box )

        secondLine = Gtk.HBox()
        thirdLine = Gtk.HBox()
        fourthLine = Gtk.HBox()
        options.pack_start(secondLine, False, False, 3 )
        label = {}
        self.paramsEntry.append( {} )
        for key in self.paramsList:
            label[key] = Gtk.Label()
            label[key].set_label(key)
            self.paramsEntry[-1][key] = Gtk.Entry()
            self.paramsEntry[-1][key].set_width_chars(5)
            secondLine.pack_start( label[key], False, False, 1)
            secondLine.pack_start( self.paramsEntry[-1][key], True, True, 1)
        if len( self.paramsEntry ) == 1:
            self.paramsEntry[-1]['type'].set_text( 'v' )
            self.paramsEntry[-1]['index'].set_text( '2' )
        elif len( self.paramsEntry ) == 2:
            self.paramsEntry[-1]['type'].set_text( 'v' )
            self.paramsEntry[-1]['index'].set_text( '2' )
        else:
            self.paramsEntry[-1]['type'].set_text( self.paramsEntry[-2]['type'].get_text() )
            self.paramsEntry[-1]['index'].set_text( self.paramsEntry[-2]['index'].get_text() )

        runButton = Gtk.Button()
        options.pack_end(runButton, False, False, 3 )
        runButton.set_label('show')
        runButton.connect( 'clicked', self.onClickShow, len(panel)-1 )

        self.images.append( Gtk.Image() )
        box.pack_start( self.images[-1], True, True, 3 )
        return

    def showImage( self, panelIndex, figname ):
        self.images[panelIndex].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( figname ) )
        if iverbose >2: print('showing figure file', figname )
        self.show_all()

    def onClickShow( self, widget, panelIndex ):
        print( '*run() invoked for panelIndex', panelIndex )
        labelMe = self.paramsEntry[panelIndex]['type'].get_text()
        indexMe = int(self.paramsEntry[panelIndex]['index'].get_text())

        if not os.path.exists('figures/'): os.makedirs('figures/')
        outname = 'figures/f_%s%s.dat'%(labelMe,indexMe)
        figname = outname.replace('.dat', '.png')
        if iverbose >2: print( 'fetching figure file', figname )

        def run():
            return self.showImage(panelIndex, figname)
        self.threads.append( threading.Thread(target=run) )
        self.threads[-1].daemon = True
        self.threads[-1].start()
        return

if __name__ == '__main__':

    print(sys.argv)
    print("hello world!")
    print(sys.float_info)
    nameset = nameConventions('auau-nexus.fzo') #the first string is the name of the zdata or zxdata file, the second number is the estimated lower number of particles, approximately no. of MC* no. of particles each MC* no. of events
    datamg = FzoDramaFactory()
    print('python script for freeze-out surface calculation class test...')
    print('the name pattern of the fzo data files is:',nameset.name)
    print('a typical data file looks like:',nameset.makeName(3))
    print('the working directory is:',os.getcwd())

    fzosurface = 'Fake'
    pdtable = 'Default'
    particleoutputlist = 'Default'
    spectrumflowset = SpectrumFlowEvaluation(fzosurface, pdtable, particleoutputlist)
    spectrumflowset.evaluateGridD3nDyptDptDphi(nameset).evaluateOutputD3nDyptDptVnVsPt().SpectrumFlowWriteOutputFile(nameset)
    plotit = PlotOnScreen()
    plotit.plotOutputFile(nameset).showPlots()

    del datamg
    del nameset
    del plotit
    print('program halt! goodbye and happy debudding!')
    exit(0)

exit(0)
