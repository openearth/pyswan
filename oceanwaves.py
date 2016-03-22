# -*- coding: utf-8 -*-
"""
Module for handling spectral ocean wave data.
  
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

def directional_spreading(dirs,pdir,ms,units='deg'):

    """
    Calculate directional spreading. 
    cdir = directional_spreading(dirs,pdir,ms,units='deg')
    where for example with pdir = 135
    dirs = np.asanyarray([0, 45, 60] + list(90+np.arange(0,19)*5.) + [210,] + list(180+np.arange(1,5)*45.))
     
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
    """
    
    def __init__(self):
        
        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.t               = []
        self.x               = []
        self.y               = []
        self.lon             = []
        self.lat             = []
        self.epsg            = [] # epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = []
        
        self.f               = [] # always implicit Hz
        self.ftype           = []
        
        self.direction       = []
        self.direction_units = '' # implicit directiontype: degrees_north or degrees_true        
        
        self.energy          = np.asarray([[[[]]]]) # [t,xy,f,dir] # last dimension 'dir' extra wrt Spec1
        self.energy_units    = ''

        #self.__dict__.update(kwargs)    
        
    def __repr__(self):

        return '<Spectrum2D ' + ' shape:[nt:' + str(self.energy.shape[0])+ \
        ',nx:'+str(self.energy.shape[1])+ ',nf:'+str(self.energy.shape[2])+ \
        ',nd:'+str(self.energy.shape[3]) + ']: "'+self.source+'">'
        
   #def Tmij(self):
   #TO DO calculate period based on various spectral moments
   
   
    
    def Hm0(self):
        """
        Integrate Hm0 (Hs) from wave spectra (using trapz rule along last dimension.)
        S = swan.Spec2()
        S.Hm0()
        """
    
        if self.energy_units[0:9] == 'm2/Hz/deg': # deg, degr, degree
            # np.abs: descending directions lead to negs
            m0 = np.trapz(np.abs(np.trapz(self.energy,self.direction)),self.f, axis=-1) # int along last dimension
        else:
            m0 = None
            print('unknown units:"',self.energy_units,'"')
        
        return 4*np.sqrt(m0)
        
    #@static
    def plot(self,fname,it=0,ix=0):        
        """
        """ 
        
        import matplotlib.pyplot as plt
        fig=plt.figure()
        fig.set_figwidth(10)
        ax = plt.axes([0.15,.15,0.8,0.8])
        ax.set_xlabel('direction')
        ax.set_ylabel('f')
        plt.pcolormesh(self.direction,self.f,self.energy[it,ix,:,:])        
        ax.set_title('[' + str(self.lat) + ','+ str(self.lon) + ']')
        plt.savefig(fname, fontsize=7, dpi=400)
        plt.close()             

class Spec1():
    """
    Class for 1D spectra. Array shape is temporal series for a list of points.
        Example
        Sp1 = Spec1()
        Sp1.f = np.linspace(0.03,.3,100)
    """
    
    def __init__(self):
        
        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.t               = []
        self.x               = []
        self.y               = []
        self.lon             = []
        self.lat             = []
        self.epsg            = [] # epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = []
        
        self.f               = [] # always Hz
        self.ftype           = []
        
        self.direction       = np.asarray([[[]]]) # [t,xy,f]
        self.direction_units = '' # implicit directiontype: degrees_north or degrees_true
        
        self.energy          = np.asarray([[[]]]) # [t,xy,f] # last dimension 'f' extra wrt Spec1
        self.energy_units    = ''
        
        self.spreading       = np.asarray([[[]]]) # [t,xy,f]
        self.spreading_units = ''
        
        #self.__dict__.update(kwargs) 
        
    def __repr__(self):
        
        return '<Spectrum1D ' + ' shape:[nt:' + str(self.energy.shape[0])+ ',nx:'+str(self.energy.shape[1])+ ',nf:'+str(self.energy.shape[2]) + ']: "'+self.source+'">'
        
   #def Tmij(self):
   #TO DO calculate period based on various spectral moments
   
    def from_jonswap(self,Hm0,Tp,
        g         = 9.81,
        gamma     = 3.3, 
        method    = 'Yamaguchi', 
        normalize = True,
        sa        = 0.07,
        sb        = 0.09):
        """
        Generate JONSWAP spectrum
        S = swan.Spec1(Hm0,Tp)
        
        """
        
        f = self.f
        
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
            
        self.energy = E
        self.energy_units = 'm2/Hz'  

        return self.energy


    def Hm0(self):
        """
        Integrate Hm0 (Hs) from wave spectra (using trapz rule along last dimension.)
        S = swan.Spec1()
        S.Hm0()
        """
    
        if self.energy_units == 'm2/Hz':
            m0 = np.trapz(self.energy,self.f, axis=-1) # int along last dimension
        else:
            m0 = None
            print('unknown units',self.energy_units)
        
        return 4*np.sqrt(m0)

    #@static
    def plot(self,fname,it=0,ix=0):        
        """
        """ 
        
        import matplotlib.pyplot as plt
        fig=plt.figure()
        fig.set_figwidth(10)
        ax = plt.axes([0.15,.15,0.8,0.8])
        plt.plot(self.f,self.energy[it,ix,:])
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel('E [' + self.energy_units + ']')
        ax.set_title('[' + str(self.lat) + ','+ str(self.lon) + ']')
        
        ax.set_xlim(0,np.max(self.f))
        plt.savefig(fname, fontsize=7, dpi=400)
        plt.close()          

        
class Spec0():
    """
    Class for spectral parameters
    """

    def __init__(self):
        
        self.buoy            = []
        self.sensor          = []
        self.source          = ''
        self.t               = []
        self.x               = []
        self.y               = []
        self.lon             = []
        self.lat             = []
        self.epsg            = [] # epsg for (x,y) <-> (lon,lat)
        self.text            = []
        self.version         = []
        
        self.Hmax            = np.asarray([[]])
        self.Hs              = np.asarray([[]])
        self.Tp              = np.asarray([[]])
        self.Tm01            = np.asarray([[]])
        self.Tm02            = np.asarray([[]])

        #self.__dict__.update(kwargs)     

        
    def __repr__(self):
        
        return '<Spectrum0D ' + ' shape:[nt:' + str(self.Hs.shape[0])+ \
        ',nx:'+str(self.Hs.shape[1])+ ']: "'+self.source+'">'
        
    def from_Spec(self,Spec):
    
        self.buoy   = Spec.buoy   
        self.sensor = Spec.sensor 
        self.source = Spec.source 
        self.t      = Spec.t      
        self.x      = Spec.x      
        self.y      = Spec.y      
        self.lon    = Spec.lon    
        self.lat    = Spec.lat    
        self.epsg   = Spec.epsg   
        self.text   = Spec.text   
        
        self.Hs     = Spec.Hm0()
    
        return self
