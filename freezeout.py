# freeze-out implementation in python
# latest update 20191015.2152
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
testmode = True

class nameConventions:

    def __init__(self,  collision_name='auau-nexus.fzo', nMinEventIndex=1, nMaxEventIndex=None ):
        self.name = collision_name
        self.eventCounter = 0
        self.nMinEventIndex=nMinEventIndex
        if nMaxEventIndex==None:
            self.nMaxEventIndex=100
        else:
            self.nMaxEventIndex=nMaxEventIndex
        self.hbarc = 1.97
        self.c = 1.0

    def makeName( self, i = 1 ):
        return self.name.replace('.fzo', '-%d.data'%i )

    def makeFZOFilename( self, i = 1 ):
        return self.name.replace('.fzo', '-%d.fzo.data'%i )

    def outputFileName( self, label1, label2='', label3='' ):
        labelme = str(label1)
        if label2 != '':
            labelme += '_'+str(label2)
            if label3 != '':
                labelme += '_'+str(label3)
        return 'zzout_'+str(labelme)+'.dat'

    @staticmethod
    def phif(px,py): #return a value between 0 to 2pi
        zero = 0.0
        two = 2.0
        three = 3.0
        if px == zero:
            if py > zero:
                phif=pi/two
            elif py < zero:
                phif=pi/two*three
            elif py == zero:
                print('fatal error in phif(), px=py=0, please check your data file')
                #print('program halt')
                #exit(1)
                phif=None
            else:
                print('fatal error in phif(), impossible option')
                print('progrm halt')
                exit(1)
        elif px < zero:
            if py > zero:
                phif=atan2(py,px) +pi
            elif py < zero:
                phif=atan2(py,px) +pi
            elif py == zero:
                phif= pi
            else:
                print('fatal error in phif(), impossible option')
                print('progrm halt')
                exit(1)
        elif px > zero:
            if py > zero:
                phif=atan2(py,px)
            elif py < zero:
                phif=atan2(py,px) +two*pi
            elif py == zero:
                phif = 0.0
            else:
                print('fatal error in phif(), impossible option')
                print('progrm halt')
                exit(1)
        return phif

class FzoDramaFactory:

    def __init__( self ):
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

    def ReadFZOsrfc( self, nameset, iEvent ):
        if iverbose > 7: print ('Subroutine FzoDramaFactory.ReadFZOsrfc() invoked!')
        outFZOfrc = FreeZeOutSurface(opt='FZOfile')
        for i in range(iEvent,iEvent+1):
            ifname = nameset.makeFZOFilename(i+1)
            if not os.path.isfile(ifname):
                print('...skipping file:',ifname)
                self.fileIsMissing = True
                return None
            self.fileIsMissing = False
            data = open( ifname, 'r' ).readlines()
            if iverbose > 8: print('reading file:',ifname)
            isFirstLine = True
            indexp = 0
            for line in data:
                linesplit = " ".join(line.split()).split(' ')
                if isFirstLine:
                    indexp += 1
                    if iverbose > 5: print(linesplit)
                    try:
                        numberOfFZOParticles = int(linesplit[0])
                        nparb_act = int(linesplit[1])
                    except ValueError:
                        print("Fatal error: the particle number does not look like an integer!")
                        exit(1)
                    numberOfItems =  len(linesplit)
                    if numberOfItems != 2:
                        print('Fatal error: the number of items of FZO info description does not match the standard. program halt!')
                        print('numberOfItems:',numberOfItems)
                        print('troublesome content:',linesplit)
                        exit(1)
                    if iverbose > 7: print('FZO particle number of current data file:',numberOfFZOParticles)
                    for iFZOParticle in range(numberOfFZOParticles):
                        memFZOsph = FreeZeOutSPH()
                        outFZOfrc.addFZOSPH(memFZOsph)
                    isFirstLine = False
                else:
                    indexp += 1
                    if iverbose > 5: print('line index:',indexp,'file batch index/total batches:',indexp,'/5')
                    numberOfItems =  len(linesplit)
                    if iverbose > 8: print('number of items in the present batch is:', numberOfItems)
                    if indexp == 2:#the batch contains numberOfFZOParticles*nparb_act+2*nparb_act items
                        if iverbose > 5: print ('indexp: ',indexp,' numberOfItems:',numberOfItems,' vs. expected value',numberOfFZOParticles*nparb_act+2*nparb_act)
                    if indexp == 3:#the batch contains 6*4(-vectors)*numberOfFZOParticles items
                        if iverbose > 5: print ('indexp: ',indexp,' numberOfItems:',numberOfItems,' vs. expected value',6*4*numberOfFZOParticles)
                        for iFZOParticle in range(numberOfFZOParticles):
                            outFZOfrc.FZOList[iFZOParticle].tauf       = float(linesplit[iFZOParticle*(4*6)+0])
                            outFZOfrc.FZOList[iFZOParticle].xf         = float(linesplit[iFZOParticle*(4*6)+1])
                            outFZOfrc.FZOList[iFZOParticle].yf         = float(linesplit[iFZOParticle*(4*6)+2])
                            outFZOfrc.FZOList[iFZOParticle].etaf       = float(linesplit[iFZOParticle*(4*6)+3])
                            outFZOfrc.FZOList[iFZOParticle].uTau       = float(linesplit[iFZOParticle*(4*6)+(4*2)+0])
                            outFZOfrc.FZOList[iFZOParticle].uX         = float(linesplit[iFZOParticle*(4*6)+(4*2)+1])
                            outFZOfrc.FZOList[iFZOParticle].uY         = float(linesplit[iFZOParticle*(4*6)+(4*2)+2])
                            outFZOfrc.FZOList[iFZOParticle].uEta       = float(linesplit[iFZOParticle*(4*6)+(4*2)+3])
                            outFZOfrc.FZOList[iFZOParticle].d3SigmaTau = float(linesplit[iFZOParticle*(4*6)+(4*3)+0])
                            outFZOfrc.FZOList[iFZOParticle].d3SigmaX   = float(linesplit[iFZOParticle*(4*6)+(4*3)+1])
                            outFZOfrc.FZOList[iFZOParticle].d3SigmaY   = float(linesplit[iFZOParticle*(4*6)+(4*3)+2])
                            outFZOfrc.FZOList[iFZOParticle].d3SigmaEta = float(linesplit[iFZOParticle*(4*6)+(4*3)+3])
                    if indexp == 4:#the batch contains 2*5*numberOfFZOParticles items
                        if iverbose > 5: print ('indexp: ',indexp,' numberOfItems:',numberOfItems,' vs. expected value',2*5*numberOfFZOParticles)
                    if indexp == 5:#the batch contains (15)*numberOfFZOParticles items
                        if iverbose > 5: print ('indexp: ',indexp,' numberOfItems:',numberOfItems,' vs. expected value',15*numberOfFZOParticles)
                        for iFZOParticle in range(numberOfFZOParticles):
                            outFZOfrc.FZOList[iFZOParticle].ef         = None
                            outFZOfrc.FZOList[iFZOParticle].pf         = None
                            outFZOfrc.FZOList[iFZOParticle].Tf         = float(linesplit[iFZOParticle*(15)+10])
                            outFZOfrc.FZOList[iFZOParticle].muBf       = float(linesplit[iFZOParticle*(15)+0])
                            outFZOfrc.FZOList[iFZOParticle].muSf       = float(linesplit[iFZOParticle*(15)+1])
        NotConsistentSPH = 0
        for iFZOParticle in range(numberOfFZOParticles-1,-1,-1):# exception check, delete must work in reversed order
            NotConsistentFlag = False
            if outFZOfrc.FZOList[iFZOParticle].Tf == 0.0 :
                NotConsistentFlag = True #remove SPH particles with vanishing freeze-out temperature
            elif outFZOfrc.FZOList[iFZOParticle].tauf == 0.0:
                NotConsistentFlag = True #remove SPH particles with vanishing initial time
            else:
                utau  = outFZOfrc.FZOList[iFZOParticle].uTau
                ux    = outFZOfrc.FZOList[iFZOParticle].uX
                uy    = outFZOfrc.FZOList[iFZOParticle].uY
                ueta  = outFZOfrc.FZOList[iFZOParticle].uEta
                uunit = utau**2-ux**2-uy**2-ueta**2/outFZOfrc.FZOList[iFZOParticle].tauf**2 #four velocity u_\mu not normalized
                if abs(uunit-1.0) > 0.001: NotConsistentFlag = True
            if NotConsistentFlag == True:
                NotConsistentSPH += 1
                outFZOfrc.deleteFZOSPHwithIndex(iFZOParticle)
        print(NotConsistentSPH,'Non-consistent FZO SPH particles have been removed!')
        return outFZOfrc


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
            self.uEta=sqrt(self.uTau**2-self.uX**2-self.uY**2-1.0)*self.tauf
            self.ef=None
            self.pf=None
            self.Tf=0.15 #freeze-out temperature scenario
            self.muBf=0.0
            self.muSf=0.0
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
        if testmode==True or opt=='Fake':
            print('randomly creating fake freeze-out surface...')
            for ifzsph in range(5):
                fzsph = FreeZeOutSPH( randomized=True )
                self.addFZOSPH(fzsph)
            print('done! ',self.numFZOSPH,' freeze-out SPH particles have been created!')

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
    def deleteFZOSPHwithIndex(self, indexSPH=None):
        if isinstance(indexSPH,int) and indexSPH<=self.numFZOSPH:
            self.FZOList.remove(self.FZOList[indexSPH])
            self.numFZOSPH -= 1
        else:
            print ('WARNING: Subroutine FreeZeOutSurface.deleteFZOSPHwithIndex() failed and ignored!')
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

#    def __init__( self, FZOsrfc='Fake', DTDMFactory='Inactive', PDTable='Default', POutList='Default', nEtaBin=None, etaMax=None, etaMin=None, nYBin=None, yMax=None, yMin=None, nPtBin=None, ptMax=None, ptMin=None, nPhiBin=None, nHarmonics=None): #python 2.7
    def __init__(self, FZOsrfc='Fake', DTDMFactory='Inactive', PDTable='Default', POutList='Default', **kwargs):
        print ('SpectrumFlowEvaluation Class instantiated')
        self.doneHeader = False
        self.gridsHaveBeenDefined = False
        self.eigenOutputDefined = False
        if testmode ==True or FZOsrfc=='Fake':
            self.FZOsrfc_dummy = True
        else:
            self.FZOsrfc_dummy = False
            self.DtDrmFctry = DTDMFactory
            if not isinstance(self.DtDrmFctry, FzoDramaFactory):
                print ('Fatal error (1) in Class SpectrumFlowEvaluation initialization, program halt!')
                exit(1)
        if testmode==True or PDTable == 'Default':
            self.PDTable = ParticleDataTable('Default')
        else:
            pass
        if testmode==True or POutList == 'Default':
            self.POutList = ParticleOutputList('Default')
        else:
            pass
        if not (isinstance(self.PDTable, ParticleDataTable) and isinstance(self.POutList, ParticleOutputList)):
            print ('Fatal error (2) in Class SpectrumFlowEvaluation initialization, program halt!')
            exit(1)
        if self.gridsHaveBeenDefined == False: #first time runner
            saved_args = locals()
            self.nEtaBin = 10 #number of grids for output, must be a factor of nEtaIntBin
            self.nEtaIntBin = 121 #integral grids must be an odd number (see https://en.wikipedia.org/wiki/Simpson%27s_rule)            self.etaMax = 6.0
            self.nYIntBin = 121
            self.etaMax = 6.0
            self.etaMin = -self.etaMax
            self.nYBin = 10
            self.nYIntBin = 121
            self.yMax = 6.0
            self.yMin = -self.yMax
            self.ptMax = 5.0
            self.ptMin = 0.0
            self.nPtBin = 10
            self.nPtIntBin = 121
            self.phiMax = 2.0*np.pi
            self.phiMin = 0.0
            self.nPhiBin = 10
            self.nPhiIntBin = 121
            self.nHarmonics = 4
            self.nPrincipleComponent = 2
            for key, value in kwargs.items():
                if   key == 'nEtaBin'    :
                    self.nEtaBin = int(value)
                    print('nEtaBin: ',self.nEtaBin )
                elif key == 'nEtaIntBin' :
                    self.nEtaIntBin = int(value)
                    print('nEtaIntBin: ',self.nEtaIntBin )
                elif key == 'etaMax'     :
                    self.etaMax = float(value)
                    print('etaMax: ',self.etaMax )
                    self.etaMin = -self.etaMax
                    print('etaMin: ',self.etaMin )
                elif key == 'etaMin'     :
                    self.etaMin = float(value)
                    print('etaMin: ',self.etaMin )
                elif key == 'nYBin'      :
                    self.nYBin = int(value)
                    print('nYBin: ',self.nYBin )
                elif key == 'nYIntBin'   :
                    self.nYIntBin = int(value)
                    print('nYIntBin: ',self.nYIntBin )
                elif key == 'yMax'       :
                    self.yMax = float(value)
                    print('yMax: ',self.yMax )
                elif key == 'yMin'       :
                    self.yMin = float(value)
                    print('yMin: ',self.yMin )
                elif key == 'nPtBin'     :
                    self.nPtBin = int(value)
                    print('nPtBin: ',self.nPtBin )
                elif key == 'nPtIntBin'  :
                    self.nPtIntBin = int(value)
                    print('nPtIntBin: ',self.nPtIntBin )
                elif key == 'ptMax'      :
                    self.ptMax = float(value)
                    print('ptMax: ',self.ptMax )
                elif key == 'ptMin'      :
                    self.ptMin = float(value)
                    print('ptMin: ',self.ptMin )
                elif key == 'nPhiBin'    :
                    self.nPhiBin = int(value)
                    print('nPhiBin: ',self.nPhiBin )
                elif key == 'nPhiIntBin' :
                    self.nPhiIntBin = int(value)
                    print('nPhiIntBin: ',self.nPhiIntBin )
                elif key == 'nHarmonics' :
                    self.nHarmonics = int(value)
                    print('nHarmonics: ',self.nHarmonics )
                elif key == 'nPrincipleComponent' :
                    self.nPrincipleComponent = int(value)
                    print('nPrincipleComponent: ',self.nPrincipleComponent )
                else:
                    print('Fatal error: keyword', key,' does not exist, program halt!')
                    exit(1)
                print ("%s == %s" %(key, value))
            if not((self.nPtIntBin % 2)==1 and (self.nPhiIntBin % 2)==1 and (self.nYIntBin % 2)==1 and (self.nEtaIntBin % 2)==1):
                print ('Fatal error: self.nXXIntBin must be odd for Simpson integral, program halt!')
                exit(1)
            if testmode == True:
                self.nPhiBin=4
                self.nPhiIntBin=21
                self.nEtaBin=4
                self.nEtaIntBin=21
                self.nYBin=4
                self.nYIntBin=21
                self.nPtBin=4
                self.nPtIntBin=21
            self.deltaEta = (self.etaMax-self.etaMin)/self.nEtaIntBin
            self.deltaY = (self.yMax-self.yMin)/self.nYIntBin
            self.deltaPt = (self.ptMax-self.ptMin)/self.nPtIntBin
            self.deltaPhi = (self.phiMax-self.phiMin)/self.nPhiIntBin
            #spectrum and flow EP
            self.D3nDyDptDphi = zeros([self.POutList.numPOutputList,self.nPtIntBin, self.nYIntBin, self.nPhiIntBin],dtype=np.float32) #np.float64 is double precision
            self.D3nDetaDptDphi = zeros([self.POutList.numPOutputList,self.nPtIntBin, self.nYIntBin, self.nPhiIntBin],dtype=np.float32)
            self.VnYPtEP = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            self.VnYEtaEP = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nEtaBin],dtype=np.float32)
            self.VnYPtCMLT = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            self.VnYEtaCMLT = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin],dtype=np.float32)
            #flow PCA
            self.Qnp = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin+1, self.nEtaBin+1],dtype=np.complex_)
            self.VnPtDelta = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin, self.nPtBin],dtype=np.complex_)
            self.VnEtaDelta = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nEtaBin, self.nEtaBin],dtype=np.complex_)

            self.gridsHaveBeenDefined = True
            print('grids definition:')
            print('nEtaBin: ',self.nEtaBin,'nEtaIntBin: ',self.nEtaIntBin,' etaMax: ',self.etaMax, ' etaMin: ',self.etaMin)
            print('nYBin: ',self.nYBin,'nYIntBin: ',self.nYIntBin,' yMax: ',self.yMax,' yMin: ',self.yMin)
            print('nPtBin: ',self.nPtBin,'nPtIntBin: ',self.nPtIntBin,' ptMax: ',self.ptMax,' ptMin: ',self.ptMin)
            print('nPhiBin: ',self.nPhiBin,'nPhiIntBin: ',self.nPhiIntBin,' phiMax: ',self.phiMax,' phiMin: ',self.phiMin)
            print('nHarmonics: ',self.nHarmonics)

    def EventAverage( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.EventAverage() invoked')
        self.VnYPtEP /= float(nameset.eventCounter)
        self.VnYEtaEP /= float(nameset.eventCounter)
        self.VnYPtCMLT /= float(nameset.eventCounter)
        self.VnYEtaCMLT /= float(nameset.eventCounter)
        return self

    def EventLoopSpectrumFlowAnalysis( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.EventLoopSpectrumFlowAnalysis() invoked')
        for iEvent in range(nameset.nMinEventIndex, nameset.nMaxEventIndex):
            print(iEvent,'/',(nameset.nMaxEventIndex-nameset.nMinEventIndex),'possible events...')
            if self.FZOsrfc_dummy == True:
                self.FZOsrfc = FreeZeOutSurface('Fake')
            else:
                self.FZOsrfc = self.DtDrmFctry.ReadFZOsrfc( nameset, iEvent )
                if self.FZOsrfc is None:
                    print('Event ',iEvent,' is empty and thus skipped!')
                    continue
            if not isinstance(self.FZOsrfc, FreeZeOutSurface):
                print ('Fatal error in Subroutine SpectrumFlowEvaluation.EventLoopSpectrumFlowAnalysis() initialization, program halt!')
                exit(1)
            print('the number of Freeze-out SPH particles is:',self.FZOsrfc.numFZOSPH)
            nameset.eventCounter += 1
            if nameset.eventCounter == 1:
                print('for the 0-th particle:')
                print('tauf,xf,yf,etaf=',self.FZOsrfc.FZOList[0].tauf,self.FZOsrfc.FZOList[0].xf,self.FZOsrfc.FZOList[0].yf,self.FZOsrfc.FZOList[0].etaf)
                print('Tf,muBf,muSf=',self.FZOsrfc.FZOList[0].Tf,self.FZOsrfc.FZOList[0].muBf,self.FZOsrfc.FZOList[0].muSf)
                print('for the 1-st particle:')
                print('tauf,xf,yf,etaf=',self.FZOsrfc.FZOList[1].tauf,self.FZOsrfc.FZOList[1].xf,self.FZOsrfc.FZOList[1].yf,self.FZOsrfc.FZOList[1].etaf)
                print('Tf,muBf,muSf=',self.FZOsrfc.FZOList[1].Tf,self.FZOsrfc.FZOList[1].muBf,self.FZOsrfc.FZOList[1].muSf)
                print('for the last particle:')
                print('tauf,xf,yf,etaf=',self.FZOsrfc.FZOList[-1].tauf,self.FZOsrfc.FZOList[-1].xf,self.FZOsrfc.FZOList[-1].yf,self.FZOsrfc.FZOList[-1].etaf)
                print('Tf,muBf,muSf=',self.FZOsrfc.FZOList[-1].Tf,self.FZOsrfc.FZOList[-1].muBf,self.FZOsrfc.FZOList[-1].muSf)

            self.evaluateGridD3nDyDptDphi(nameset).evaluateGridD3nDetaDptDphi(nameset)
            self.evaluateOutputD3nDyptDptVnVsPt(EPCalculated=False).evaluateOutputD3nDyptDptVnVsPt(EPCalculated=True)
            self.buildVarianceMatrix(nameset)

        print('The total number of events included in the present analysis is: ',nameset.eventCounter)
        self.normalizeVarianceMatrix( nameset ).eigenValueSolver(nPrincipleComponent=2).pcaFlowWriteOutputFile(nameset)
        self.EventAverage( nameset )
        return self

    def header( self, nameset ):
        f = open( nameset.outputFileName('vn-pt'), 'w' )
        f.write('pt'+' '+'vn'+'\n')
        f.close()
        f = open( nameset.outputFileName('vn-eta'), 'w' )
        f.write('eta'+' '+'vn'+'\n')
        f.close()
        f = open( nameset.outputFileName('pca-pt'), 'w' )
        f.write('lambdaPC '+' '.join([' '.join(['{:f}'.format(lambdaPC) for lambdaPC in givenHarmonics]) for givenHarmonics in self.eigenPtValue[:,self.POutList.numPOutputList,:]])+'\n')
        if iverbose > 7: print('pt '+' '.join([' '.join(['n{0!s}p{1!s}'.format(givenHarmonics,pComponent) for pComponent in range(self.nPrincipleComponent)]) for givenHarmonics in range(self.nHarmonics) ])+'\n')
        f.write('pt '+' '.join([' '.join(['n{0!s}pc{1!s}'.format(givenHarmonics,pComponent) for pComponent in range(self.nPrincipleComponent)]) for givenHarmonics in range(self.nHarmonics) ])+'\n')
        f.close()
        f = open( nameset.outputFileName('pca-eta'), 'w' )
        f.write('lambdaPC '+' '.join([' '.join(['{:f}'.format(lambdaPC) for lambdaPC in givenHarmonics]) for givenHarmonics in self.eigenEtaValue[:,self.POutList.numPOutputList,:]])+'\n')
        f.write('eta '+' '.join([' '.join(['n{0!s}pc{1!s}'.format(givenHarmonics,pComponent) for pComponent in range(self.nPrincipleComponent)]) for givenHarmonics in range(self.nHarmonics) ])+'\n')
        f.close()
        self.doneHeader = True

    def SpectrumFlowWriteOutputFile( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.SpectrumFlowOutput() invoked')
        if not self.doneHeader:
            self.header( nameset )
        f = open( nameset.outputFileName('vn-pt'), 'a' )
        for iPtBin in range(self.nPtBin):
            pt = self.ptMin+(iPtBin+0.5)*(self.ptMax-self.ptMin)
            f.write(str(pt)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYPtEP[:, self.POutList.numPOutputList, iPtBin]])+'\n')
            if iverbose > 8: print(str(pt)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYPtEP[:, self.POutList.numPOutputList, iPtBin]])+'\n')
        f.close()
        f = open( nameset.outputFileName('vn-eta'), 'a' )
        for iEtaBin in range(self.nEtaBin):
            eta = self.etaMin+(iEtaBin+0.5)*(self.etaMax-self.etaMin)
            f.write(str(eta)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYEtaEP[:, self.POutList.numPOutputList, iEtaBin]])+'\n')
            if iverbose > 8: print(str(eta)+' '+' '.join(['{:f}'.format(givenHarmonics) for givenHarmonics in self.VnYEtaEP[:, self.POutList.numPOutputList, iEtaBin]])+'\n')
        f.close()
        print('The obtained results on Vn(EP) have been saved to the file(s).')
        return self

    def pcaFlowWriteOutputFile( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.pcaFlowWriteOutputFile() invoked')
        if not self.doneHeader:
            self.header( nameset )
        f = open( nameset.outputFileName('pca-pt'), 'a' )
        for j in range(self.nPtBin):
            pt = self.ptMin+(j+0.5)*(self.ptMax-self.ptMin)
            f.write(str(pt)+' '+' '.join([' '.join(['{:f}'.format(pComponent) for pComponent in givenHarmonics]) for givenHarmonics in self.eigenPtVector[:,self.POutList.numPOutputList,:,j]])+'\n') # https://stackoverflow.com/questions/17870612/printing-a-two-dimensional-array-in-python and https://pyformat.info/
            if iverbose > 8: print(str(pt)+' '+' '.join([' '.join(['{:f}'.format(pComponent) for pComponent in givenHarmonics]) for givenHarmonics in self.eigenPtVector[:,self.POutList.numPOutputList,:,j]])+'\n')
        f.close()
        f = open( nameset.outputFileName('pca-eta'), 'a' )
        for j in range(self.nEtaBin):
            eta = self.etaMin+(j+0.5)*(self.etaMax-self.etaMin)
            f.write(str(eta)+' '+' '.join([' '.join(['{:f}'.format(pComponent) for pComponent in givenHarmonics]) for givenHarmonics in self.eigenEtaVector[:,self.POutList.numPOutputList,:,j]])+'\n')
        f.close()
        print('The obtained results on Vn(PCA) have been saved to the file(s).')
        return self

    def buildVarianceMatrix( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.buildVarianceMatrix() invoked')
        QnpTmpa = zeros([self.nHarmonics, self.POutList.numPOutputList+1, self.nPtBin+1, self.nEtaBin+1],dtype=np.complex_)
        if ( ((self.nPtIntBin-1)%self.nPtBin) == 0 and ((self.nPtIntBin-1)/self.nPtBin)%2 == 0):
            print('Fatal error in PCA analysis')
            print('condition not satisfied: self.nPtIntBin - 1 = even number * self.nPtBin')
            print('self.nPtIntBin= ',self.nPtIntBin,' self.nPtBin= ',self.nPtBin)
            print('((self.nPtIntBin-1)%self.nPtBin)= ',((self.nPtIntBin-1)%self.nPtBin))
            print('((self.nPtIntBin-1)/self.nPtBin)%2= ',((self.nPtIntBin-1)/self.nPtBin)%2)
            print('program halt!')
            exit(1)
        if ( ((self.nEtaIntBin-1)%self.nEtaBin) == 0 and ((self.nEtaIntBin-1)/self.nEtaBin)%2 == 0):
            print('Fatal error in PCA analysis')
            print('condition not satisfied: self.nEtaIntBin - 1 = even number * self.nEtaBin')
            print('self.nEtaIntBin= ',self.nEtaIntBin,' self.nEtaBin= ',self.nEtaBin)
            print('((self.nEtaIntBin-1)%self.nEtaBin)= ',((self.nEtaIntBin-1)%self.nEtaBin))
            print('((self.nEtaIntBin-1)/self.nEtaBin)%2= ',((self.nEtaIntBin-1)/self.nEtaBin)%2)
            print('program halt!')
            exit(1)
        for iHarmonics in range(self.nHarmonics):
            for iParticleSpeciesOutput in range(self.POutList.numPOutputList):
                for iPtBin in range(self.nPtBin):
                    iPtIntBinMin=iPtBin*(self.nPtIntBin-1)/self.nPtBin
                    iPtIntBinMax=iPtIntBinMin+(self.nPtIntBin-1)/self.nPtBin+1
                    if iverbose > 8: print('iPtIntBinMin,iPtIntBinMax= ',iPtIntBinMin,iPtIntBinMax)
                    for iEtaBin in range(self.nEtaBin):
                        iEtaIntBinMin=iEtaBin*(self.nEtaIntBin-1)/self.nEtaBin
                        iEtaIntBinMax=iEtaIntBinMin+(self.nEtaIntBin-1)/self.nEtaBin+1
                        if iverbose > 8: print('iEtaIntBinMin,iEtaIntBinMax= ',iEtaIntBinMin,iEtaIntBinMax,'... ',end='')
                        sum0 = 0j
                        for iPtIntBin in range(iPtIntBinMin,iPtIntBinMax):
                            iPtIntBinEff=iPtIntBin-iPtIntBinMin
                            nPtIntBinEff=iPtIntBinMax-iPtIntBinMin
                            if (iPtIntBinEff == 0 or iPtIntBinEff == nPtIntBinEff-1):
                                fator0 = 1.0
                            elif iPtIntBinEff%2 == 0:
                                fator0 = 4.0
                            elif iPtIntBinEff%2 == 1:
                                fator0 = 2.0
                            sum1 = 0j
                            for iEtaIntBin in range(iEtaIntBinMin,iEtaIntBinMax):
                                if iverbose > 8 and iPtIntBin == iPtIntBinMin: print(iEtaIntBin,end='')
                                iEtaIntBinEff=iEtaIntBin-iEtaIntBinMin
                                nEtaIntBinEff=iEtaIntBinMax-iEtaIntBinMin
                                if (iEtaIntBinEff == 0 or iEtaIntBinEff == nEtaIntBinEff-1):
                                    fator1 = 1.0
                                elif iEtaIntBinEff%2 == 0:
                                    fator1 = 4.0
                                elif iEtaIntBinEff%2 == 1:
                                    fator1 = 2.0
                                sum2 = 0j
                                for iPhiIntBin in range(self.nPhiIntBin):
                                    fator2=3.0 #periodic function
                                    rphi=self.phiMin+(iPhiIntBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiIntBin)
                                    tempa = exp(1j*iHarmonics*rphi)*self.D3nDetaDptDphi[iParticleSpeciesOutput,iPtIntBin,iEtaIntBin,iPhiIntBin]
                                    tempa *= fator2
                                    sum2 += tempa*fator2
                                sum1 += sum2*fator1
                            sum0 += sum1*fator0
                        if iverbose > 8:print('')
                        sum0 *= (self.deltaPt/3.0)*(self.deltaEta/3.0)*(self.deltaPhi/3.0)
                        IntervalPt=(self.ptMax-self.ptMin)/self.nPtBin
                        IntervalEta=(self.etaMax-self.etaMin)/self.nEtaBin
                        QnpTmpa[iHarmonics,iParticleSpeciesOutput,iPtBin,iEtaBin] = sum0/(2.0*pi*IntervalPt*IntervalEta)
                        QnpTmpa[iHarmonics,iParticleSpeciesOutput,iPtBin,self.nPtBin] = sum0/(2.0*pi*IntervalPt)
                        QnpTmpa[iHarmonics,iParticleSpeciesOutput,self.nEtaBin,iEtaBin] = sum0/(2.0*pi*IntervalEta)
                        QnpTmpa[iHarmonics,self.POutList.numPOutputList,iPtBin,iEtaBin] += sum0/(2.0*pi*IntervalPt*IntervalEta)
                        QnpTmpa[iHarmonics,self.POutList.numPOutputList,iPtBin,self.nPtBin] += sum0/(2.0*pi*IntervalPt)
                        QnpTmpa[iHarmonics,self.POutList.numPOutputList,self.nEtaBin,iEtaBin] += sum0/(2.0*pi*IntervalEta)

        QnpTmpPtDiagnal = diag( QnpTmpa[0,self.POutList.numPOutputList,:self.nPtBin,self.nEtaBin] )
        if iverbose>8: print('QnpTmpa[0,self.POutList.numPOutputList,:self.nPtBin,self.nEtaBin]:\n',QnpTmpa[0,self.POutList.numPOutputList,:self.nPtBin,self.nEtaBin])
        if iverbose>8: print('QnpTmpPtDiagnal:\n',QnpTmpPtDiagnal)
        self.VnPtDelta += -1.0*QnpTmpPtDiagnal/(2.0*pi*IntervalPt) # broadcast, the word here is broadcast
        if iverbose>8: print('self.VnPtDelta:\n',self.VnPtDelta)
        self.VnPtDelta += QnpTmpa[:,:,:self.nPtBin,None,self.nEtaBin] * conjugate(QnpTmpa[:,:,None,:self.nPtBin,self.nEtaBin])
        if iverbose>8: print('self.VnPtDelta:\n',self.VnPtDelta)

        QnpTmpEtaDiagnal = diag( QnpTmpa[0,self.POutList.numPOutputList,self.nPtBin,:self.nEtaBin] )
        self.VnEtaDelta += -1.0*QnpTmpEtaDiagnal/(2.0*pi*IntervalEta) # broadcast, the word here is broadcast
        self.VnEtaDelta += QnpTmpa[:,:,self.nPtBin,:self.nEtaBin,None] * conjugate(QnpTmpa[:,:,self.nPtBin,None,:self.nEtaBin])
        self.Qnp += QnpTmpa
        return self

    def normalizeVarianceMatrix(self, nameset):
        if iverbose > 2: print ('Subroutine SpectrumFlowEvaluation.normalizeVarianceMatrix() invoked!')
        print('To normalize the variance matrix based on a total of',nameset.eventCounter,'events has been scanned.')
        self.Qnp /= nameset.eventCounter # normalization for event number
        self.VnPtDelta /= nameset.eventCounter
        self.VnEtaDelta /= nameset.eventCounter
        if iverbose > 4:
            print('deltaEta',self.deltaEta,'deltaPt',self.deltaPt)
            for iHarmonics in range(self.nHarmonics):
                print('Qnp[',iHarmonics,']:')
                print(self.Qnp[Harmonics,self.POutList.numPOutputList])
        if iverbose > 4:
            for iHarmonics in range(self.nHarmonics):
                print('<QnpQnp> VnPtDelta[',iHarmonics,']:')
                print(self.VnPtDelta[Harmonics,self.POutList.numPOutputList])
                print('<QnpQnp> VnEtaDelta[',iHarmonics,']:')
                print(self.VnEtaDelta[Harmonics,self.POutList.numPOutputList])
        self.VnPtDelta += -1.0*self.Qnp[:,:,:self.nPtBin,None,self.nEtaBin] * conjugate(self.Qnp[:,:,None,:self.nPtBin,self.nEtaBin])
        self.VnEtaDelta += -1.0*self.Qnp[:,:,self.nPtBin,:self.nEtaBin,None] * conjugate(self.Qnp[:,:,self.nPtBin,None,:self.nEtaBin])
        return self

    def eigenValueSolver(self, nPrincipleComponent=3):
        if iverbose > 2: print ('Subroutine SpectrumFlowEvaluation.eigenValueSolver() invoked!')
        if self.eigenOutputDefined == False:
            self.nPrincipleComponent = nPrincipleComponent
            self.eigenPtValue = zeros([self.nHarmonics,self.POutList.numPOutputList+1,self.nPrincipleComponent],dtype=np.complex_)
            self.eigenPtVector = zeros([self.nHarmonics,self.POutList.numPOutputList+1,self.nPrincipleComponent,self.nPtBin],dtype=np.complex_)
            self.eigenEtaValue = zeros([self.nHarmonics,self.POutList.numPOutputList+1,self.nPrincipleComponent],dtype=np.complex_)
            self.eigenEtaVector = zeros([self.nHarmonics,self.POutList.numPOutputList+1,self.nPrincipleComponent,self.nEtaBin],dtype=np.complex_)
            self.eigenOutputDefined = True
        for iHarmonics in range(self.nHarmonics):
            eigenPtValueTemp, eigenPtVectorTemp = linalg.eig(self.VnPtDelta[iHarmonics,self.POutList.numPOutputList])
            eigenPtPairTemp = [(eigenPtValueTemp[i], eigenPtVectorTemp[:,i]) for i in range(len(eigenPtValueTemp))]
            eigenPtPairTemp.sort(reverse = True, key = lambda elem : abs(real(elem[0])))
            for i in range(self.nPrincipleComponent):
                self.eigenPtValue[iHarmonics,self.POutList.numPOutputList,i] = eigenPtPairTemp[i][0] # damn! list and array are different https://www.pythoncentral.io/the-difference-between-a-list-and-an-array/
                for j in range(self.nPtBin):
                    self.eigenPtVector[iHarmonics,self.POutList.numPOutputList,i,j] = eigenPtPairTemp[i][1][j]
            for i in range(self.nPrincipleComponent):
                for j in range(self.nPtBin):
                    jpt = j+1
                    self.eigenPtVector[iHarmonics,self.POutList.numPOutputList,i,j] *= sqrt(abs(self.eigenPtValue[iHarmonics,self.POutList.numPOutputList,i]))/self.Qnp[0,self.POutList.numPOutputList,jpt,self.nEtaBin]
            eigenEtaValueTemp, eigenEtaVectorTemp = linalg.eig(self.VnEtaDelta[iHarmonics,self.POutList.numPOutputList])
            eigenEtaPairTemp = [(eigenEtaValueTemp[i], eigenEtaVectorTemp[:,i]) for i in range(len(eigenEtaValueTemp))]
            eigenEtaPairTemp.sort(reverse = True, key = lambda elem : abs(real(elem[0])))
            for i in range(self.nPrincipleComponent):
                self.eigenEtaValue[iHarmonics,self.POutList.numPOutputList,i] = eigenEtaPairTemp[i][0]
                for j in range(self.nEtaBin):
                    self.eigenEtaVector[iHarmonics,self.POutList.numPOutputList,i,j] = eigenEtaPairTemp[i][1][j]
            for i in range(self.nPrincipleComponent):
                for j in range(self.nEtaBin):
                    jeta = j+1
                    self.eigenEtaVector[iHarmonics,self.POutList.numPOutputList,i,j] *= sqrt(abs(self.eigenEtaValue[iHarmonics,self.POutList.numPOutputList,i]))/self.Qnp[0,self.POutList.numPOutputList,self.nPtBin,jeta]
            if iverbose > 4:
                print('eigenPtValue[',iHarmonics,'], eigenPtVector[',iHarmonics,']:')
                print(self.eigenPtValue[iHarmonics,self.POutList.numPOutputList], self.eigenPtVector[iHarmonics,self.POutList.numPOutputList])
                print('eigenPtValue[',iHarmonics,'], eigenPtVector[',iHarmonics,']:')
                print(self.eigenEtaValue[iHarmonics,self.POutList.numPOutputList], self.eigenEtaVector[iHarmonics,self.POutList.numPOutputList])
        return self



    def evaluateGridD3nDyDptDphi( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.evaluateGridD3nDyDptDphi() invoked')
        print(self.FZOsrfc.numFZOSPH,end="")
        for iescolTmp in range(self.FZOsrfc.numFZOSPH):
            print('.',end="")
            for iParticleSpeciesOutput in range(self.POutList.numPOutputList):
                for iPtIntBin in range(self.nPtIntBin):
                    rptTmp=self.ptMin+(iPtIntBin+0.5)*(self.ptMax-self.ptMin)/float(self.nPtIntBin)
                    for iYIntBin in range(self.nYIntBin):
                        ryTmp=self.yMin+(iYIntBin+0.5)*(self.yMax-self.yMin)/float(self.nYIntBin)
                        for iPhiIntBin in range(self.nPhiIntBin):
                            rphiTmp=self.phiMin+(iPhiIntBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiIntBin)
                            tmpDis=self.fD3nDptDyDphi(nameset,rpt=rptTmp,ry=ryTmp,rphi=rphiTmp,iescol=iescolTmp,ipid=iParticleSpeciesOutput)
                            self.D3nDyDptDphi[iParticleSpeciesOutput,iPtIntBin,iYIntBin,iPhiIntBin]+=tmpDis
        print(' ')
        return self

    def evaluateGridD3nDetaDptDphi( self, nameset ):
        print ('Subroutine SpectrumFlowEvaluation.evaluateGridD3nDetaDptDphi() invoked')
        print(self.FZOsrfc.numFZOSPH,end="")
        for iescolTmp in range(self.FZOsrfc.numFZOSPH):
            print('.',end="")
            for iParticleSpeciesOutput in range(self.POutList.numPOutputList):
                for iPtIntBin in range(self.nPtIntBin):
                    rptTmp=self.ptMin+(iPtIntBin+0.5)*(self.ptMax-self.ptMin)/float(self.nPtIntBin)
                    for iEtaIntBin in range(self.nEtaIntBin):
                        retaTmp=self.etaMin+(iEtaIntBin+0.5)*(self.etaMax-self.etaMin)/float(self.nEtaIntBin)
                        for iPhiIntBin in range(self.nPhiIntBin):
                            rphiTmp=self.phiMin+(iPhiIntBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiIntBin)
                            tmpDis=self.fD3nDptDetaDphi(nameset,rpt=rptTmp,reta=retaTmp,rphi=rphiTmp,iescol=iescolTmp,ipid=iParticleSpeciesOutput)
                            self.D3nDetaDptDphi[iParticleSpeciesOutput,iPtIntBin,iEtaIntBin,iPhiIntBin]+=tmpDis
        print(' ')
        return self

    def getiPtIntBin(self,rpt): #the output index is from 0 to nPtBin-1, -1 is reserved for out-of-range input
        iPtIntBin=int(float(self.nPtIntBin)*(rpt-self.ptMin)/(self.ptMax-self.ptMin))
        if iPtIntBin < 0: #pt cannot be negative
            print('Fatal error in getipt()')
            print('rpt:',rpt,'ptMax:',self.ptMax,'ptMin:',self.ptMin,'nPtIntBin:',self.nPtIntBin)
            print('program halt')
            exit(1)
        if iPtIntBin >= self.nPtIntBin: iPtIntBin=-1 #if out of the range then one should let it go
        return iPtIntBin

    def getiEtaIntBin(self,reta): #the output index is from 0 to nPtBin-1, -1 is reserved for out-of-range input
        iEtaIntBin=int(float(self.nEtaIntBin)*(reta-self.etaMin)/(self.etaMax-self.etaMin))
        if iEtaIntBin < 0: iEtaIntBin=-1
        if iEtaIntBin >= self.nEtaIntBin: iEtaIntBin=-1 #if out of the range then one should let it go
        return iEtaIntBin

    def evaluateOutputD3nDyptDptVnVsPt( self, EPCalculated=False ):
        print ('Subroutine evaluateOutputD3nDyptDptVnVsPt() invoked')
# 0 - pt, 0.5 (0a5) - particle species, 1 - phi, 2 - rapidity (y)
        sum0 = 0.0
        if EPCalculated==False: #event plane not calculated, so evaluate them!
            self.phievt = zeros([self.nHarmonics],dtype=np.float32)
            self.psib   = zeros([self.nHarmonics],dtype=np.float32)
            self.psia   = zeros([self.nHarmonics],dtype=np.float32)
            self.psiptb = zeros([self.nHarmonics],dtype=np.float32)
            self.psipta = zeros([self.nHarmonics],dtype=np.float32)
        fase    = zeros([self.nHarmonics],dtype=np.float32)
        psisin  = zeros([self.nHarmonics],dtype=np.float32)
        psicos  = zeros([self.nHarmonics],dtype=np.float32)
        psibsin = zeros([self.nHarmonics],dtype=np.float32)
        psibcos = zeros([self.nHarmonics],dtype=np.float32)
        psiasin = zeros([self.nHarmonics],dtype=np.float32)
        psiacos = zeros([self.nHarmonics],dtype=np.float32)

        for iPtBin in range(self.nPtBin):
            rpt=self.ptMin+(iPtBin+0.5)*(self.ptMax-self.ptMin)/float(self.nPtBin)
            iPtIntBin=self.getiPtIntBin(rpt)
            if (iPtBin == 0 or iPtBin == self.nPtBin-1):
                fator0 = 1.0
            elif iPtBin%2 == 0:
                fator0 = 4.0
            elif iPtBin%2 == 1:
                fator0 = 2.0
            sum0a5 = 0.0
            for iHarmonics in range(self.nHarmonics):
                soma0a5cos = zeros([self.nHarmonics],dtype=np.float32)
            for iParticleSpeciesOutput in range(self.POutList.numPOutputList):
                sum1 = 0.0
                for iHarmonics in range(self.nHarmonics):
                    soma1sin    = zeros([self.nHarmonics],dtype=np.float32)
                    soma1cos    = zeros([self.nHarmonics],dtype=np.float32)
                    soma1bsin   = zeros([self.nHarmonics],dtype=np.float32)
                    soma1bcos   = zeros([self.nHarmonics],dtype=np.float32)
                    soma1asin   = zeros([self.nHarmonics],dtype=np.float32)
                    soma1acos   = zeros([self.nHarmonics],dtype=np.float32)
                for iPhiIntBin in range(self.nPhiIntBin):
                    rphi=self.phiMin+(iPhiIntBin+0.5)*(self.phiMax-self.phiMin)/float(self.nPhiIntBin)
                    fator1=3.0 #periodic function
                    sum2 = 0.0
                    for iHarmonics in range(self.nHarmonics):
                        soma2sin    = zeros([self.nHarmonics],dtype=np.float32)
                        soma2cos    = zeros([self.nHarmonics],dtype=np.float32)
                        soma2bsin   = zeros([self.nHarmonics],dtype=np.float32)
                        soma2bcos   = zeros([self.nHarmonics],dtype=np.float32)
                        soma2asin   = zeros([self.nHarmonics],dtype=np.float32)
                        soma2acos   = zeros([self.nHarmonics],dtype=np.float32)
                    for iYIntBin in range(self.nYIntBin):
                        ry=self.yMin+(iYIntBin+0.5)*(self.yMax-self.yMin)/float(self.nYIntBin)
                        if (ry < 0.0):
                            fase[iHarmonics]   = self.psib[iHarmonics]
                        else:
                            fase[iHarmonics]   = self.psia[iHarmonics]
                        if (iYIntBin == 1 or iYIntBin == self.nYIntBin-1):
                            fator2 = 1.0
                        elif iYIntBin%2 == 0:
                            fator2 = 4.0
                        elif iYIntBin%2 == 1:
                            fator2 = 2.0
                        tmp2    = self.D3nDyDptDphi[iParticleSpeciesOutput,iPtIntBin,iYIntBin,iPhiIntBin]
#                        if tmp2 == 0.0:
#                            print('WARNING: vanishing tmp2 observed, iParticleSpeciesOutput=',iParticleSpeciesOutput,' iPtBin=',iPtBin,' iYBin=',iYBin,' iPhiBin=',iPhiBin)
#                            print('program halt!')
#                            exit(1)
                        sum2   += tmp2*fator2
                        for iHarmonics in range(self.nHarmonics):
                            soma2sin[iHarmonics]+=tmp2*fator2*sin(float(iHarmonics)*(rphi-fase[iHarmonics]))
                            soma2cos[iHarmonics]+=tmp2*fator2*cos(float(iHarmonics)*(rphi-fase[iHarmonics]))
                            if ry > 0.0:
                                soma2bsin[iHarmonics]+=tmp2*fator2*sin(float(iHarmonics)*(rphi-fase[iHarmonics]))
                                soma2bcos[iHarmonics]+=tmp2*fator2*cos(float(iHarmonics)*(rphi-fase[iHarmonics]))
                            else:
                                soma2asin[iHarmonics]+=tmp2*fator2*sin(float(iHarmonics)*(rphi-fase[iHarmonics]))
                                soma2acos[iHarmonics]+=tmp2*fator2*cos(float(iHarmonics)*(rphi-fase[iHarmonics]))
                    sum1   += sum2*fator1
                    for iHarmonics in range(self.nHarmonics):
                        soma1sin[iHarmonics]  += soma2sin[iHarmonics]*fator1
                        soma1cos[iHarmonics]  += soma2cos[iHarmonics]*fator1
                        soma1bsin[iHarmonics] += soma2bsin[iHarmonics]*fator1
                        soma1bcos[iHarmonics] += soma2bcos[iHarmonics]*fator1
                        soma1asin[iHarmonics] += soma2asin[iHarmonics]*fator1
                        soma1acos[iHarmonics] += soma2acos[iHarmonics]*fator1
                sum1   *= self.deltaY/3.0
                sum1   *= self.deltaPhi/3.0
                for iHarmonics in range(self.nHarmonics):
                    soma1sin[iHarmonics]  *= self.deltaY/3.0
                    soma1sin[iHarmonics]  *= self.deltaPhi/3.0
                    soma1cos[iHarmonics]  *= self.deltaY/3.0
                    soma1cos[iHarmonics]  *= self.deltaPhi/3.0
                    soma1bsin[iHarmonics] *= self.deltaY/3.0
                    soma1bsin[iHarmonics] *= self.deltaPhi/3.0
                    soma1bcos[iHarmonics] *= self.deltaY/3.0
                    soma1bcos[iHarmonics] *= self.deltaPhi/3.0
                    soma1asin[iHarmonics] *= self.deltaY/3.0
                    soma1asin[iHarmonics] *= self.deltaPhi/3.0
                    soma1acos[iHarmonics] *= self.deltaY/3.0
                    soma1acos[iHarmonics] *= self.deltaPhi/3.0
                sum0a5    += sum1 #only to sum up different specials, therefore no fator0
                sum0      += sum1*fator0 #this summation in pt is taken place also for every species
                for iHarmonics in range(self.nHarmonics):
                    soma0a5cos[iHarmonics] += soma1cos[iHarmonics]
                    psisin[iHarmonics]     += soma1sin[iHarmonics]*fator0
                    psicos[iHarmonics]     += soma1cos[iHarmonics]*fator0
                    psibsin[iHarmonics]    += soma1bsin[iHarmonics]*fator0
                    psibcos[iHarmonics]    += soma1bcos[iHarmonics]*fator0
                    psiasin[iHarmonics]    += soma1asin[iHarmonics]*fator0
                    psiacos[iHarmonics]    += soma1acos[iHarmonics]*fator0

                self.VnYPtEP[0,iParticleSpeciesOutput,iPtBin] += sum1 #multiplicity for given particle species
                for iHarmonics in range(1,self.nHarmonics):
                    if sum1 != 0:
                        self.VnYPtEP[iHarmonics,iParticleSpeciesOutput,iPtBin] += psicos[iHarmonics]/sum1 #vn for given particle species
                    else:
                        print('Fatal error: vanishing sum1 observed, iParticleSpeciesOutput=',iParticleSpeciesOutput,' iPtBin=',iPtBin)
                        print('program halt!')
                        exit(1)

            self.VnYPtEP[0,self.POutList.numPOutputList,iPtBin] += sum0a5 #multiplicity, all particle species
            for iHarmonics in range(1,self.nHarmonics):
                if sum0a5 != 0:
                    self.VnYPtEP[iHarmonics,self.POutList.numPOutputList,iPtBin] += psicos[iHarmonics]/sum0a5 #vn, all particle species
                else:
                    print('WARNING: vanishing sum0a5 observed, iPtBin=',iPtBin)
                    print('program halt!')
                    exit(1)

        sum0   *= self.deltaPt/3.0

        for iHarmonics in range(self.nHarmonics):
            psisin[iHarmonics]  *= self.deltaPt/3.0
            psicos[iHarmonics]  *= self.deltaPt/3.0
            psibsin[iHarmonics] *= self.deltaPt/3.0
            psibcos[iHarmonics] *= self.deltaPt/3.0
            psiasin[iHarmonics] *= self.deltaPt/3.0
            psiacos[iHarmonics] *= self.deltaPt/3.0

        if EPCalculated==False:
            for iHarmonics in range(1,self.nHarmonics):#avoid dividing by zero
                self.phievt[iHarmonics] = nameConventions.phif(psicos[iHarmonics],psisin[iHarmonics])/float(iHarmonics)
                self.psib[iHarmonics]   = nameConventions.phif(psibcos[iHarmonics],psibsin[iHarmonics])/float(iHarmonics)
                self.psia[iHarmonics]   = nameConventions.phif(psiacos[iHarmonics],psiasin[iHarmonics])/float(iHarmonics)

        return self

    def fD3nDptDyDphi(self,nameset,rpt=None,ry=None,rphi=None,iescol=None,ipid=None):
        saved_args = locals()
        if str(saved_args).count('None') != 0 :
            print('fD3nDptDyDphi() input fatal error, program halt!')
            exit(1)

        detetmin = 0.170             #the highest possible value
        Tfcut    = 0.100             #the lowest allowed value
        Tfzo = min(detetmin,self.FZOsrfc.FZOSPH(iescol).Tf)
        if Tfzo < Tfcut:
            d3ndptdydphi=0.0
            return d3ndptdydphi
        ub = self.FZOsrfc.FZOSPH(iescol).muBf
        us = self.FZOsrfc.FZOSPH(iescol).muSf
        b  = self.PDTable.EleParticle(ipid).b
        m  = self.PDTable.EleParticle(ipid).m
        g  = self.PDTable.EleParticle(ipid).g
        ssq  = self.PDTable.EleParticle(ipid).s
        sinal = self.PDTable.EleParticle(ipid).sgl
        #print('iescol,Tf',iescol,self.FZOsrfc.FZOSPH(iescol).Tf)
        #print(type(iescol),type(self.FZOsrfc.FZOSPH(iescol).Tf))
        #print('Tfzo,m,g,sinal=',Tfzo,m,g,sinal)
        #print(type(Tfzo),type(m),type(g),type(sinal))
        #print('b*ub+ssq*us',b,ub,ssq,us)
        #print(type(b),type(ub),type(ssq),type(us))
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
        if (1.0-v2) < 0.0:
            print('NaN in d3n! program halt!')
            print('(1.d0-v2)=',(1.0-v2))
            print('utau,ux,uy,ueta=',utau,ux,uy,ueta)
            print('tau,x,y,eta=',tau,x,y,eta)
            print('Tfzo=',Tfzo)
            exit(1)
        gamma  = 1.0/sqrt(1.0-v2)

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

    def fD3nDptDetaDphi(self,nameset,rpt=None,reta=None,rphi=None,iescol=None,ipid=None):
        saved_args = locals()
        if str(saved_args).count('None') != 0 :
            print('fD3nDptDetaDphi() input fatal error, program halt!')
            exit(1)

        detetmin = 0.170             #the highest possible value
        Tfcut    = 0.100             #the lowest allowed value
        Tfzo = min(detetmin,self.FZOsrfc.FZOSPH(iescol).Tf)
        if Tfzo < Tfcut:
            d3ndptdetadphi=0.0
            return d3ndptdetadphi
        ub = self.FZOsrfc.FZOSPH(iescol).muBf
        us = self.FZOsrfc.FZOSPH(iescol).muSf
        b  = self.PDTable.EleParticle(ipid).b
        m  = self.PDTable.EleParticle(ipid).m
        g  = self.PDTable.EleParticle(ipid).g
        ssq  = self.PDTable.EleParticle(ipid).s
        sinal = self.PDTable.EleParticle(ipid).sgl
        #print('iescol,Tf',iescol,self.FZOsrfc.FZOSPH(iescol).Tf)
        #print(type(iescol),type(self.FZOsrfc.FZOSPH(iescol).Tf))
        #print('Tfzo,m,g,sinal=',Tfzo,m,g,sinal)
        #print(type(Tfzo),type(m),type(g),type(sinal))
        #print('b*ub+ssq*us',b,ub,ssq,us)
        #print(type(b),type(ub),type(ssq),type(us))
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
        if (1.0-v2) < 0.0:
            print('NaN in d3n! program halt!')
            print('(1.d0-v2)=',(1.0-v2))
            print('utau,ux,uy,ueta=',utau,ux,uy,ueta)
            print('tau,x,y,eta=',tau,x,y,eta)
            print('Tfzo=',Tfzo)
            exit(1)
        gamma  = 1.0/sqrt(1.0-v2)

# 4-vetor p^{\mu} (txyz-coord); lab frame
        p0=sqrt((rpt*cosh(reta))**2+m**2)
#        p0=mt*cosh(ry)
        px=rpt*cos(rphi)
        py=rpt*sin(rphi)
        pz=rpt*sinh(reta)
#        pz=mt*sinh(ry)

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
        d3ndptdetadphi=taa.real
        return d3ndptdetadphi

    def evaluateOutputD3nDetaDpTDeta():
        print ('Subroutine evaluateOutputDn3DetapTdpT() invoked')
        return self

    def evaluateOutputDnDeta():
        print ('Subroutine evaluateOutputDnDeta() invoked')
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

    def readOutputFile( self, nameset, iHarmonics=2, nPrincipleComponent=None, iPrincipleComponent=0 ):
        if nPrincipleComponent == None:
            print('fatal error in variable nPrincipleComponent input, program halt!')
            exit(1)
        ptep = []
        vneppt = []
        etaep = []
        vnepeta = []
        data = open( nameset.outputFileName('vn-pt'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ') # reorganize the obtained line into strings separated by single blank ' '
            if linesplit[0] == 'pt':
                lambdaPt= 1
                continue
            else:
                ptep += [ float(linesplit[0]) ]
                vneppt += [ float(linesplit[iHarmonics+1]) ]
        data = open( nameset.outputFileName('vn-eta'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
            if linesplit[0] == 'eta':
                lambdaEta= 2
                continue
            else:
                etaep += [ float(linesplit[0]) ]
                vnepeta += [ float(linesplit[iHarmonics+1]) ]
        ptep = array(ptep)
        vneppt = array(vneppt)
        etaep = array(etaep)
        vnepeta = array(vnepeta)
        if iverbose > 8: print('ptep, vneppt, etaep, vnepeta\n',ptep, vneppt, etaep, vnepeta)

        ptpca = []
        rept = []
        impt = []
        etapca = []
        reeta = []
        imeta = []
        iPrincipleComponent = 0
        data = open( nameset.outputFileName('pca-pt'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ') # reorganize the obtained line into strings separated by single blank ' '
            if linesplit[0] == 'pt':
                continue
            elif linesplit[0] == 'lambdaPC':
                lambdaPt= complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])
                continue
            ptpca += [ float(linesplit[0]) ]
            rept += [ real(complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])) ]
            impt += [ imag(complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])) ]
        data = open( nameset.outputFileName('pca-eta'), 'r' ).readlines()
        for line in data:
            linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
            if linesplit[0] == 'eta':
                continue
            elif linesplit[0] == 'lambdaPC':
                lambdaEta= complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])
                continue
            etapca += [ float(linesplit[0]) ]
            reeta += [ real(complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])) ]
            imeta += [ imag(complex(linesplit[nPrincipleComponent*iHarmonics+iPrincipleComponent+1])) ]
        ptpca = array(ptpca)
        rept = array(rept)
        impt = array(impt)
        etapca = array(etapca)
        reeta = array(reeta)
        imeta = array(imeta)
        if iverbose > 8: print('ptpca, rept, impt, etapca, reeta, imeta\n',ptpca, rept, impt, etapca, reeta, imeta)

        return ptep, vneppt, ptpca, rept, etaep, vnepeta, etapca, reeta

    def plotOutputFile( self, nameset, nHarmonics=None, nPrincipleComponent=None, iPrincipleComponent=0 ):
        if nHarmonics==None or nPrincipleComponent==None:
            print('fatal error in variable nHarmonics or nPrincipleComponent, program halt!')
            exit(1)
        for iHarmonics in range(nHarmonics):
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

            x1, y1, x2, y2, x3, y3, x4, y4 = self.readOutputFile( nameset, iHarmonics, nPrincipleComponent )
            labelMe = 'n'
            indexMe = iHarmonics

            subFig1.plot(x1[:], y1[:], '-b', lw=.5, label='EP')
            subFig2.plot(x2[:], y2[:], '-b', lw=.5, label='PCA' )
            subFig3.plot(x3[:], y3[:], '-r', lw=.5, label='EP' )
            subFig4.plot(x4[:], y4[:], '-r', lw=.5, label='PCA' )

            subFig1.set_ylabel(r'$v_{:}$'.format(indexMe), fontsize=16)
            subFig2.set_ylabel(r'$v_{:}$'.format(indexMe), fontsize=16)
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
            self.paramsEntry[-1]['type'].set_text( 'n' )
            self.paramsEntry[-1]['index'].set_text( '2' )
        elif len( self.paramsEntry ) == 2:
            self.paramsEntry[-1]['type'].set_text( 'n' )
            self.paramsEntry[-1]['index'].set_text( '0' )
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
    nameset = nameConventions('auau-xell.fzo',nMaxEventIndex=4) #the first string is the name of the zdata or zxdata file, the second number is the estimated lower number of particles, approximately no. of MC* no. of particles each MC* no. of events
    datamg = FzoDramaFactory()
    print('python script for freeze-out surface calculation class test...')
    print('the name pattern of the fzo data files is:',nameset.name)
    print('a typical data file looks like:',nameset.makeName(3))
    print('the working directory is:',os.getcwd())

    fzosurface = 'FZOfile' #'FZOfile' or 'Fake'
    pdtable = 'Default'
    particleoutputlist = 'Default'
    spectrumflowset = SpectrumFlowEvaluation(fzosurface, datamg, pdtable, particleoutputlist,nEtaBin=11,nEtaIntBin=121)
    spectrumflowset.EventLoopSpectrumFlowAnalysis(nameset).SpectrumFlowWriteOutputFile(nameset)
    plotit = PlotOnScreen()
    plotit.plotOutputFile(nameset, spectrumflowset.nHarmonics, spectrumflowset.nPrincipleComponent).showPlots()

    del plotit
    del spectrumflowset
    del datamg
    del nameset
    print('program halt! goodbye and happy debudding!')
    exit(0)

exit(0)
