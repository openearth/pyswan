# -*- coding: utf-8 -*-
"""
Module for handling spectral ocean wave data. We adhere to the following convention:
 * directions:   nautical convention: CF units degrees_true
 * time:         UTC
 * coordinates:  speherical WGS84
 * energy units: [m2/Hz/deg]
 * methods for spectral parameters Hs(), Tm01, Tm02 
   use the trapezozidal rule to integrate the spectrum.
  
"""

"""
$HeadURL:  $
$Id: m $
"""
__version__ = "$Revision: $"

#  Copyright notice
#   --------------------------------------------------------------------
#   Copyright (C) 2016 TU Delft/Van Oord
#       Gerben J. de Boer, <g.j.deboer@tudelft.nl>/<gerben.deboer@vanoord.com>
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this library.  If not, see <http://www.gnu.org/licenses/>.
#   --------------------------------------------------------------------

# This tool is part of <a href="http://www.OpenEarth.eu">OpenEarthTools</a>.
# OpenEarthTools is an online collaboration to share and manage data and
# programming tools in an open source, version controlled environment.
# Sign up to recieve regular updates of this function, and to contribute
# your own tools.

import numpy as np
import datetime
debug = False

def jonswap(f,Hm0,Tp,
    g         = 9.81,
    gamma     = 3.3, 
    method    = 'Yamaguchi', 
    normalize = True,
    sa        = 0.07,
    sb        = 0.09):
    """
    Generate 1D JONSWAP spectrum in [m2/Hz']
    E = ow.jonswap(f,Hm0,Tp,**kwargs)
    By default E is nornalized to the integral of E == Hm0.
    
    """
    
    # method    = 'Yamaguchi'; # 'Goda'        

    def sigma(f,fpeak,sa,sb):
        s = np.ones(f.size)*sa
        s[f > fpeak] = sb
        return s
        
    # Pierson-Moskowitz
    if method=='Yamaguchi':
        alpha = 1/(0.06533*gamma**0.8015 + 0.13467)/16; # Yamaguchi (1984), used in SWAN
    elif method=='Goda':
        alpha = 1/(0.23+0.03*gamma-0.185*(1.9+gamma)**-1)/16; # Goda

    pm  = alpha*Hm0**2*Tp**-4*f**-5*np.exp(-1.25*(Tp*f)**-4);
    
    # apply JONSWAP shape
    E = pm*gamma**np.exp(-0.5*(Tp*f-1)**2./sigma(f,1/Tp,sa,sb)**2);
    #E(np.isnan(E))=0

    if normalize:
        corr = Hm0**2/(16*np.trapz(E,f))
        E = E*corr
        
    return E

def directional_spreading(dirs,pdir,ms,units='deg'):

    """
    Calculate directional spreading. 
    >> cdir = directional_spreading(dirs,pdir,ms,units='deg')
    where for example with pdir = 135
    >> dirs = np.asanyarray([0, 45, 60] + list(90+np.arange(0,19)*5.) + \
              [210,] + list(180+np.arange(1,5)*45.))
     
    """
    from math import gamma
    dirs = np.asarray(dirs) # also accept list

    if units=='deg':
        dirs = dirs*np.pi/180
        pdir = pdir*np.pi/180
    elif units=='rad':
        dirs = dirs + 0. # ensure real
        pass
    else:
        raise('unknown units')

    A1 = (2.**ms) * (gamma(ms/2+1))**2 / (np.pi * gamma(ms+1))
    cdir = dirs*0
    for i,v in enumerate(dirs):
        acos = np.cos((dirs[i] - pdir))
        if acos > 0:
            cdir[i] = A1*np.max(acos**ms, 1e-10)
    if units=='deg':
        cdir = cdir*np.pi/180
    elif units=='rad':
        pass
    else:
        raise('unknown units')    
        
    # check for ill-sampling
    intz = np.trapz(cdir,dirs)*180/np.pi
    if not(abs(intz-1)< 1e-6):
        print('integral not 1 for ms=',ms,':',intz) 

    return cdir
   
class Spec2():
    """
    Class for 2D spectra. Array shape is temporal series for a list of points.
    Example:
    >> Sp2 = Spec2(f=np.linspace(0.03,.3,100),direction=np.arange(0,24)*15) 
    generates 2D energy array with nans
    """
    
    def __init__(self,f=[np.nan],direction=[np.nan],**kwargs):
        
        self.f               = f # always implicit Hz
        self.direction       = direction

        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.t               = [np.nan]
        self.x               = [np.nan]
        self.y               = [np.nan]
        self.lon             = [np.nan]
        self.lat             = [np.nan]
        self.epsg            = []# epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = [] 
        
        self.energy_units    = ''
        self.direction_units = '' # implicit directiontype: degrees_north or degrees_true, radians??      

        self.__dict__.update(kwargs)    
        
        nt = len(self.t)
        nx = len(self.x)
        nf = len(self.f)
        nd = len(self.direction)
        
        self.energy          = np.asarray(np.nan*np.zeros((nt,nx,nf,nd)))  # [t,xy,f,dir] or [xy,f,dir] or [f,dir] # last dimension 'dir' extra wrt Spec1
        
        energy = np.asarray(kwargs.pop('energy',self.energy))
        if energy.shape==(len(self.t),len(self.x),len(self.f),len(self.direction)):
            self.energy = np.asarray(energy)
        else:
            raise Exception('dimensions E '+str(energy.shape)+' do not match t,x,f,direction '+str((len(self.t),len(self.x),len(self.f),len(self.direction))))        
        
    def __repr__(self):
    
        if   len(self.energy.shape)==2:
            txt = '<Spectrum2D ' + ' shape:[nt: 0' + \
            ',nx: 0 '+ ',nf:'+str(self.energy.shape[0])+ \
            ',nd:'+str(self.energy.shape[1]) + ']: "'+self.source+'">'
        elif len(self.energy.shape)==3:
            txt = '<Spectrum2D ' + ' shape:[nt: 0' + \
            ',nx:'+str(self.energy.shape[0])+ ',nf:'+str(self.energy.shape[1])+ \
            ',nd:'+str(self.energy.shape[2]) + ']: "'+self.source+'">'        
        elif len(self.energy.shape)==4:
            txt = '<Spectrum2D ' + ' shape:[nt:' + str(self.energy.shape[0])+ \
            ',nx:'+str(self.energy.shape[1])+ ',nf:'+str(self.energy.shape[2])+ \
            ',nd:'+str(self.energy.shape[3]) + ']: "'+self.source+'">'        

        return txt
        
   #def Tmij(self):
   #TO DO calculate period based on various spectral moments
   
    
    def Hm0(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec2(), S.Hm0()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            # np.abs: descending directions lead to negs
            m0 = np.trapz(np.abs(np.trapz(self.energy,self.direction)),self.f, axis=-1) # int along last dimension
        else:
            m0 = None
            print('unknown units:"',self.energy_units,'"')
        
        return 4*np.sqrt(m0)

    def Tm01(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec2(), S.Tm01()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            # np.abs: descending directions lead to negs
            m0 = np.trapz(np.abs(np.trapz(self.energy,self.direction)),self.f, axis=-1) # int along last dimension
            m1 = np.trapz(np.abs(np.trapz(self.energy,self.direction))*self.f,self.f, axis=-1) # int along last dimension
            Tm = m0/m1
        else:
            Tm = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tm


    def Tm02(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec2(), S.Tm02()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            # np.abs: descending directions lead to negs
            m0 = np.trapz(np.abs(np.trapz(self.energy,self.direction)),self.f, axis=-1) # int along last dimension
            m2 = np.trapz(np.abs(np.trapz(self.energy,self.direction))*self.f**2,self.f, axis=-1) # int along last dimension
            Tm = np.sqrt(m0/m2)
        else:
            Tm = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tm

    def Tp(self):
        """
        Get peak period: S = swan.Spec2(), S.Tp()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            Tp=1/self.f[np.argmax(np.max(self.energy,axis=-1))]
        else:
            Tp = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tp    

    def pdir(self):
        """
        Get peak period: S = swan.Spec2(), S.pdir()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            pdir = self.direction[np.argmax(np.max(self.energy,axis=-2))]
        else:
            pdir = None
            print('unknown units:"',self.energy_units,'"')
        
        return pdir         
        
    def from_jonswap(self,Hm0,Tp,pdir,ms,**kwargs):
        """
        Generate 2D JONSWAP spectrum
        >> Sp2 = swan.Spec2(f=f,t=t,x=x,direction=d)
        >> Sp2.from_jonswap(Hm0,Tp)
        
        """

        energy1 = jonswap(self.f,Hm0,Tp,**kwargs)
        #self.energy = np.tile(Sp1.energy,[len(self.direction),1]).T
        for it in range(len(self.t)):
            for ix in range(len(self.x)):
                for iff,v in enumerate(self.f):
                    cdir = directional_spreading(self.direction,pdir,ms) # ms[iff]
                    self.energy[it,ix,iff,:] = cdir*energy1[iff];        
        
        self.energy_units = 'm2/Hz/deg'

        return self
        
    #@static
    def plot(self,fname=None,it=0,ix=0):        
        """
        plot 2D spectrum (to file)
        """ 
        
        f = np.tile(self.f        ,[len(self.direction),1]).T
        d = np.tile(self.direction,[len(self.f),1])
        fx,lx = f*np.cos(90-d*np.pi/180),'f'
        fy,ly = f*np.sin(90-d*np.pi/180),'f'
        
        #fx,lx = self.direction,'direction'
        #fy,ly = self.f        ,'f'
        
        import matplotlib.pyplot as plt
        if not fname is None:
            fig=plt.figure()
            fig.set_figwidth(10)
        ax = plt.axes([0.15,.15,0.8,0.8])
        ax.set_xlabel(lx)
        ax.set_ylabel(ly)
        #if   len(self.energy.shape)==2:        
        #    plt.pcolormesh(fx,fy,self.energy[:,:])        
        #elif   len(self.energy.shape)==3:        
        #    plt.pcolormesh(fx,fy,self.energy[ix,:,:])        
        #elif   len(self.energy.shape)==4:        
        plt.pcolormesh(fx,fy,self.energy[it,ix,:,:])        
        ax.set_title(str(self.t[it])+'  [' + str(self.lat[ix]) + ' °N,'+ str(self.lon[ix]) + ' °E]')
        plt.axis('equal')
        
        if not fname is None:
            plt.savefig(fname, fontsize=7, dpi=400)
            plt.close()
        else:
            return ax

class Spec1():
    """
    Class for 1D spectra. Array shape is temporal series for a list of points.
    Example:
    >> Sp1 = Spec1(f=np.linspace(0.03,.3,100)) 
    generates 1D energy array with nans
    """
    
    def __init__(self,f=[np.nan],**kwargs):
        
        self.f               = f  # always implicit Hz
        
        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.t               = [np.nan]
        self.x               = [np.nan]
        self.y               = [np.nan]
        self.lon             = [np.nan]
        self.lat             = [np.nan]
        self.epsg            = []# epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = []         
        
        self.energy_units    = ''
        self.direction_units = '' # implicit directiontype: degrees_north or degrees_true
        self.spreading_units = ''
        
        self.__dict__.update(kwargs) # get f for making E,direction,spreading
        
        nt = len(self.t)
        nx = len(self.x)
        nf = len(self.f)
        
        self.energy          = np.asarray(np.nan*np.zeros((nt,nx,nf)))  # [t,xy,f] # last dimension 'f' extra wrt Spec1
        self.direction       = np.asarray(np.nan*np.zeros((nt,nx,nf)))  # [t,xy,f] or [xy,f] or [f]
        self.spreading       = np.asarray(np.nan*np.zeros((nt,nx,nf)))  # [t,xy,f]
        
        energy = np.asarray(kwargs.pop('energy',self.energy))
        if energy.shape==(len(self.t),len(self.x),len(self.f)):
            self.energy = np.asarray(energy)
        else:
            raise Exception('dimensions E '+str(energy.shape)+' do not match t,x,f '+str((len(self.t),len(self.x),len(self.f))))
            
        direction = np.asarray(kwargs.pop('direction',self.direction))
        if direction.shape==(len(self.t),len(self.x),len(self.f)):
            self.direction = np.asarray(direction)
        else:
            raise Exception('dimensions E '+str(direction.shape)+' do not match t,x,f '+str((len(self.t),len(self.x),len(self.f))))

        spreading = np.asarray(kwargs.pop('spreading',self.spreading))
        if spreading.shape==(len(self.t),len(self.x),len(self.f)):
            self.spreading = np.asarray(spreading)
        else:
            raise Exception('dimensions E '+str(spreading.shape)+' do not match t,x,f '+str((len(self.t),len(self.x),len(self.f))))            
        
    def __repr__(self):
        
        if   len(self.energy.shape)==1:
            txt = '<Spectrum1D ' + ' shape:[nt: 0' + ',nx: 0'+ ',nf:'+str(self.energy.shape[0]) + ']: "'+self.source+'">'
        elif   len(self.energy.shape)==2:
            txt = '<Spectrum1D ' + ' shape:[nt: 0' + ',nx:'+str(self.energy.shape[0])+ ',nf:'+str(self.energy.shape[1]) + ']: "'+self.source+'">'
        elif len(self.energy.shape)==3:
            txt = '<Spectrum1D ' + ' shape:[nt:' + str(self.energy.shape[0])+ ',nx:'+str(self.energy.shape[1])+ ',nf:'+str(self.energy.shape[2]) + ']: "'+self.source+'">'        
        
        return txt
        
        
   #def Tmij(self):
   #TO DO calculate period based on various spectral moments
   
    def from_jonswap(self,Hm0,Tp,**kwargs):
        """
        Generate 1D JONSWAP spectrum
        Sp1 = swan.Spec1(f=f,t=t,x=x)
        Sp1.from_jonswap(Hm0,Tp)
        
        """
        
        for it in range(len(self.t)):
            for ix in range(len(self.x)):
                self.energy[it,ix] = jonswap(self.f,Hm0,Tp,**kwargs)
        
        self.energy_units = 'm2/Hz'  

        return self

    def Hm0(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec1(), S.Hm0()
        """
    
        if self.energy_units == 'm2/Hz':
            m0 = np.trapz(self.energy,self.f, axis=-1) # int along last dimension
        else:
            m0 = None
            print('unknown units',self.energy_units)
        
        return 4*np.sqrt(m0)


    def Tm01(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec1(), S.Tm01()
        """
    
        if self.energy_units[0:9] == 'm2/Hz':
            # np.abs: descending directions lead to negs
            m0 = np.trapz(self.energy,self.f       , axis=-1) # int along last dimension
            m1 = np.trapz(self.energy*self.f,self.f, axis=-1) # int along last dimension
            Tm = m0/m1
        else:
            Tm = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tm


    def Tm02(self):
        """
        Integrate Hm0 (Hs) from wave spectra: S = swan.Spec1(), S.Tm02()
        """
    
        if self.energy_units[0:9] == 'm2/Hz':
            # np.abs: descending directions lead to negs
            m0 = np.trapz(self.energy,self.f          , axis=-1) # int along last dimension
            m2 = np.trapz(self.energy*self.f**2,self.f, axis=-1) # int along last dimension
            Tm = np.sqrt(m0/m2)
        else:
            Tm = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tm


    def Tp(self):
        """
        Get peak period: S = swan.Spec1(), S.Tp()
        """
    
        if self.energy_units[0:9] == 'm2/Hz':
            Tp = 1/self.f[np.argmax(self.energy)]
        else:
            Tp = None
            print('unknown units:"',self.energy_units,'"')
        
        return Tp        

    #@static
    def plot(self,fname=None,it=0,ix=0):        
        """
        plot 1D spectrum (to file)
        """ 
        
        import matplotlib.pyplot as plt
        if not fname is None:
            fig=plt.figure()
            fig.set_figwidth(10)
        ax = plt.axes([0.15,.15,0.8,0.8])
        #if   len(self.energy.shape)==1:
        #    plt.plot(self.f,self.energy[:])
        #elif   len(self.energy.shape)==2:
        #    plt.plot(self.f,self.energy[ix,:])
        #elif   len(self.energy.shape)==3:
        plt.plot(self.f,self.energy[it,ix,:])
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel('E [' + self.energy_units + ']')
        ax.set_title(str(self.t[it])+'  [' + str(self.lat[ix]) + ' °N,'+ str(self.lon[ix]) + ' °E]')
        
        if not fname is None:
            plt.savefig(fname, fontsize=7, dpi=400)
            plt.close()
        else:
            return ax

        
class Spec0():
    """
    Class for spectral parameters
    """

    def __init__(self,**kwargs):
        
        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.x               = []
        self.y               = []
        self.lon             = []
        self.lat             = []
        self.epsg            = [] # epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = []
        
        self.__dict__.update(kwargs)           
        
        self.t               = np.asarray([[]])
        self.Hs              = np.asarray([[]])
        self.Tp              = np.asarray([[]])
        self.Tm01            = np.asarray([[]])
        self.Tm02            = np.asarray([[]])
        self.pdir            = np.asarray([[]])
        self.ms              = np.asarray([[]])
        
    def __repr__(self):
    
        if type(self.Hs) is type(1.) or type(self.Hs) is type(1):
        
            return '<Spectrum0D  Hs='+str(self.Hs)+\
            ' Tp='+str(self.Tp)+\
            ' @ '+str(self.t)+\
            ' : "'+self.source+'">'            
        
        else:
        
            return '<Spectrum0D ' + ' shape:[nt,nx:' + str(self.Hs.shape)+']: "'+self.source+'">'

        
    def from_Spec(self,Spec):
    
        self.buoy   = Spec.buoy   
        self.sensor = Spec.sensor 
        self.source = Spec.source 
        self.x      = Spec.x      
        self.y      = Spec.y      
        self.lon    = Spec.lon    
        self.lat    = Spec.lat    
        self.epsg   = Spec.epsg   
        self.text   = Spec.text   
        
        self.t      = Spec.t      
        self.Hs     = Spec.Hm0()
        self.Tp     = Spec.Tp()
        self.Tm01   = Spec.Tm01()
        self.Tm02   = Spec.Tm02()
        self.pdir   = Spec.pdir()
       #self.ms     = Spec.pdir() # TO DO fit jonswap with Tp
    
        return self
        
        
    #@static
    def plot(self,fname=None,it=0,ix=0):        
        """
        plot parameter time series (to file)
        """ 
        
        import matplotlib.pyplot as plt
        if not fname is None:
            fig=plt.figure()
            fig.set_figwidth(10)
            
        if type(self.Hs) is type(1.) or type(self.Hs) is type(1):
        
            print('scalar: not plotted.')
            ax = None        
            
        else:
            
            ax = []
            ax.append(plt.axes([0.1,0.75,0.8,0.2], axisbg='w'))
            plt.plot_date(self.t,self.Hs,'-')
            ax[-1].set_ylabel('Hs [m]')

            ax.append(plt.axes([0.1,0.5,0.8,0.2], axisbg='w'))
            plt.plot_date(self.t,self.Tp,'-')
            ax[-1].set_ylabel('Tp [s]')

            ax.append(plt.axes([0.1,0.25,0.8,0.2], axisbg='w'))
            plt.plot_date(self.t,self.pdir,'-')
            ax[-1].set_ylabel('pdir [deg]')

            ax.append(plt.axes([0.1,0.0,0.8,0.2], axisbg='w'))
            plt.plot_date(self.t,self.ms,'-')
            ax[-1].set_ylabel('ms [-]')

        
        if not fname is None:
            plt.savefig(fname, fontsize=7, dpi=400)
            plt.close()
        else:
            return ax     
