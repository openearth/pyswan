# -*- coding: utf-8 -*-
"""
Module for handling SWAN input/output, using oceanwwaves package:
>> import oceanwaves as ow
  
"""

__version__ = "$Revision: WIP$" # https://git-scm.com/book/en/v2/Customizing-Git-Git-Attributes

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

# This tool is part of <a href="https://github.com/openearth/">OpenEarth</a>.
# OpenEarth is an online collaboration to share and manage data and
# programming tools in an open source, version controlled environment.
# Sign up to recieve regular updates of this function, and to contribute
# your own tools.

import numpy as np
import datetime
import oceanwaves as ow
import os.path

debug = False
   
#@static
def from_swan_header(f, source='from_swan'):
    """
    """
    
    self = {}
    
    self["source"] = source
    self["text"] = []
    self["version"] = f.readline().split()[1].strip()
    raw = f.readline()
    while raw.strip()[0] == '$':
        self["text"].append(raw.split('$')[1].rstrip('\n'))
        raw = f.readline()
    
    # when

    if raw.split()[0].strip().upper() == 'TIME':
        timecoding = int(f.readline().split()[0].strip())
        raw = f.readline()
    else:
        timecoding = None

    # where

    xy  = raw.split()[0].strip()
    nxy = int(f.readline().split()[0])
    self["x"]   = []
    self["y"]   = []
    self["lon"] = []
    self["lat"] = []
    
    for i in range(nxy):
        raw = f.readline()
        if debug:
            print('POINT ', i,': ',raw.rstrip())
        if xy.upper() == 'LONLAT':
            self["lon"].append(float(raw.split()[0]))
            self["lat"].append(float(raw.split()[1]))
            self["x"].append(None)
            self["y"].append(None)
        else:
            self["x"].append(float(raw.split()[0]))
            self["y"].append(float(raw.split()[1]))
            self["lon"].append(None)
            self["lat"].append(None)

    # frequency

    self["f"] = []
    self["ftype"] = f.readline().split()[0].strip()
    
    if self["ftype"].upper() == 'CDIR':
        self["ftype"] = 'cartesian'
    else: # NDIR
        self["ftype"] = 'nautical'
        
    nf = int(f.readline().split()[0])
    for i in range(nf):
        raw = f.readline()
        self["f"].append(float(raw.split()[0]))
    self["f"] = np.asarray(self["f"])
        
    return self, timecoding

def from_file2D(f,source='from_swan'):
    """
    Load 2D spectral series file from SWAN spectral file/buffer.
    
    Matrix 'energy' has 4 dimensions:
    [time x locations x frequencies x directions]
    
    f = fopen(r'c:\temp\myfile.spc')
    S2 = sean.from_swan(f,source='myfile.spc') # optionally remember source
    f.close()
    """           
    self = ow.Spec2()
    
    D, timecoding = from_swan_header(f)
    for k in D.keys():
        self.__setattr__(k, D[k])    
    self.source = source
    
    # direction

    self.direction = []
    direction_type = f.readline().split()[0].strip()
    nf = int(f.readline().split()[0])
    for i in range(nf):
        raw = f.readline()
        self.direction.append(float(raw.split()[0]))
    self.direction = np.asarray(self.direction)
    
    if direction_type.upper() == 'CDIR':
        self.direction_units = 'degrees_true'  # CF convention
    else:    
        self.direction_units = 'degrees_north' # CF convention
    
    # what
    
    f.readline() # QUANT
    number_of_quantities = int(f.readline().split()[0]) # 1
    quantity_names = []
    quantity_units = []
    quantity_exception_values = []
    for j in range(number_of_quantities):
        quantity_names.append(f.readline().split()[0].strip())
        quantity_units.append(f.readline().split()[0].strip())
        quantity_exception_values.append(float(f.readline().split()[0].strip()))
    self.energy_units    = quantity_units[0]
    
    # data
    
    raw = 'ini'
    En4,t=[],[]
    raw = f.readline()
    while not raw == '':
        if timecoding == None: # time is optional
            t.append(np.nan)
        else:
            tstr = raw.split()[0].strip()
            if timecoding==1:
                t.append(datetime.datetime.strptime(tstr,'%Y%m%d.%H%M%S'))
            else:
                t.append(tstr)
                print('timecoding unknown: ',timecoding)
            if debug:
                print('TIME',len(t),': ',t[-1])
            raw = f.readline()
        En3 = []
        for i in range(len(self.x)):
            loc = raw.split()[0]
            if loc =='FACTOR':
                factor = float(f.readline().split()[0].strip())
                if debug:
                    print('   POINT ',i,' FACTOR:',factor)
            else:
                print('ERROR, only FACTOR not implemented not ZERO or NODATA ',raw)
            raw = f.readline()
            En2 = []
            for j in range(len(self.f)):
                En1 = []
                for k in range(len(self.direction)):
                    En1.append(float(raw.split()[k])*factor)
                En2.append(En1)
                raw = f.readline()
                if raw == '':
                    break
            En3.append(En2)
        En4.append(En3)
    En4         = np.asarray(En4)
    self.t      = np.asarray(t)
    self.energy = np.ma.masked_array(En4,En4==quantity_exception_values[0])

    return self 
    
#@static
def to_file2D(self,f,timecoding=1):
    """
    Save 2D spectral series file to SWAN spectral file/buffer.
    
    S2 = ow.Spec2() 
    # fill S2
    f = fopen(r'c:\temp\newfile.s2d')
    S2.to_file2D(f)
    f.close()      
    Use attribute timecoding=None to optionally write a stationary file
    i.e. without time. When time is nan, also stationary file will be written.
    
    """      
    
    f.writelines('SWAN ' + str(self.version) + '\n')
    it = 0
    for it,t in enumerate(self.text):
        f.writelines('$' + t + '\n')
        
    # TIME    
    if not timecoding is None:
        if type(self.t[0]) is type(datetime.datetime(1997,1,1)):
            f.writelines('TIME\n')
            #f.writelines(str(len(self.x)) + '\n') # can SWAN handle this?, that would be very convenient for pre-allocation
            f.writelines(str(min(1,len(self.x))) + '\n')    
            
    # LOC
    if self.lon[0] == None or np.isnan(self.lon[0]):
        f.writelines('LOCATIONS\n')
        f.writelines(str(len(self.x)) + '\n')
        for i in range(len(self.x)):
            f.write("{:g} {:g}\n".format(self.x[i],self.y[i]))
        nxy = len(self.x)
    else:
        f.writelines('LONLAT\n')
        f.writelines(str(len(self.lon)) + '\n')
        for i in range(len(self.lon)):
            f.write("{:f} {:f}\n".format(self.lon[i],self.lat[i]))
        nxy = len(self.lon)
        
    # FREQ
    try:
        if self.ftype.lower() == 'absolute':
            f.writelines('AFREQ\n')
        else:    
            f.writelines('RFREQ\n')
    except:
            f.writelines('AFREQ\n')
    f.writelines(str(len(self.f))+ '\n')
    for freq in self.f:
        f.writelines("{:g}\n".format(freq))

    # DIR
    if self.direction_units.lower() == 'degrees_true':
        f.writelines('CDIR\n')
    elif self.direction_units.lower() == 'degrees_north':
        f.writelines('NDIR\n')
    else:    
        f.writelines('NDIR\n')
        print('no units defined, assumed NDIR')
    f.writelines(str(len(self.direction))+ '\n')
    for dir in self.direction:
        f.writelines("{:g}\n".format(dir))
        
    # QUANT
    f.writelines('QUANT \n')
    f.writelines('1                                       number of quantities in table\n')
    
    f.writelines('VaDens                                  variance densities in m2/Hz/degr\n')
    f.writelines(self.energy_units + '                              unit\n')
    f.writelines('NaN                                     exception value' + '\n')
    
    for it,t in enumerate(self.t):
        if not timecoding is None:
            if type(self.t[0]) is type(datetime.datetime(1997,1,1)):
                if debug:
                    print('time ',self.t[it])        
                f.writelines("{:%Y%m%d.%H%M%S}                         date and time\n".format(self.t[it]))
        for ix in range(nxy):
            f.writelines('FACTOR \n1\n'.format(ix+1)) # LOCATION
            for ifr,fr in enumerate(self.f):
                for idr,dr in enumerate(self.direction):
                    f.writelines('{:g} '.format(self.energy[it,ix,ifr,idr]))
                f.writelines('\n')    
                    
    return                
    
#@static
def from_file1D(f,source='from_swan'):
    """
    Load 1D spectral series file from SWAN spectral file/buffer.
    
    Matrices 'energy','direction','spreading' have 3 dimensions:
    [time x locations x frequencies]
    
    f = fopen(r'c:\temp\myfile.s1d')
    S1 = swan.from_file1D(f,source='myfile.s1d') # optionally remember source
    f.close()
    """   

    self = ow.Spec1()
    
    D, timecoding = from_swan_header(f)
    for k in D.keys():
        self.__setattr__(k, D[k])
    self.source = source
    
    # what

    f.readline() # QUANT
    number_of_quantities = int(f.readline().split()[0]) # 3
    quantity_names = []
    quantity_units = []
    quantity_exception_values = []
    for j in range(number_of_quantities):
        quantity_names.append(f.readline().split()[0].strip())
        quantity_units.append(f.readline().split()[0].strip())
        quantity_exception_values.append(float(f.readline().split()[0].strip()))
    self.energy_units    = quantity_units[0]
    self.direction_units = quantity_units[1]
    self.spreading_units = quantity_units[2]
    if quantity_names[1].upper() == 'CDIR':
        self.direction_units = 'degrees_true'  # CF convention
    else:    
        self.direction_units = 'degrees_north' # CF convention
    if quantity_names[2] == 'DSPRD' or quantity_names[2] == 'DEGR':
        self.spreading_units = 'degrees' # DSPRD or DEGR
    else:    
        self.spreading_units = '1'       # DSPRP or POWER         
        
    # data

    raw = 'ini'
    En3,Th3,Sp3,t = [],[],[],[]
    if debug:
        print('timecoding',timecoding)
        print('TIMES')
    
    raw = f.readline()
    if debug:
        print('-1-',raw)
    while not raw == '':
        if timecoding == None: # time is optional
            t.append(np.nan)
        else:
            tstr = raw.split()[0].strip()
            if timecoding==1:
                t.append(datetime.datetime.strptime(tstr,'%Y%m%d.%H%M%S'))
            else:
                t.append(tstr)
                print('timecoding unknown: ',timecoding)
            if debug:
                print(len(t),': ',t[-1])
            raw = f.readline()
            if debug:
                print('-2-',raw)
        En2,Th2,Sp2 = [],[],[]
        for i in range(len(self.x)):
            raw = f.readline()
            if debug:
                print('-3-',raw)
#  TO DO: can be NODATA
            En1,Th1,Sp1 = [],[],[]
            for j in range(len(self.f)):
                En1.append(float(raw.split()[0]))
                Th1.append(float(raw.split()[1]))
                Sp1.append(float(raw.split()[2]))
                raw = f.readline()
                if raw == '':
                    break
            En2.append(En1)
            Th2.append(Th1)
            Sp2.append(Sp1)
        En3.append(En2)
        Th3.append(Th2)
        Sp3.append(Sp2)
    En3 = np.asarray(En3)
    Th3 = np.asarray(Th3)
    Sp3 = np.asarray(Sp3)
    
    self.t         = np.asarray(t)
    self.energy    = np.ma.masked_array(En3,En3==quantity_exception_values[0])
    self.direction = np.ma.masked_array(Th3,Th3==quantity_exception_values[1])
    self.spreading = np.ma.masked_array(Sp3,Sp3==quantity_exception_values[2])
   
    return self
    
#@static
def to_file1D(self,f,timecoding=1):
    """
    Save 1D spectral series file to SWAN spectral file/buffer.
    
    S1 = ow.Spec1() 
    # fill S1
    f = fopen(r'c:\temp\newfile.s1d')
    S1.to_file1D(f)
    f.close()      
    Use attribute timecoding=None to optionally write a stationary file
    i.e. without time. When time is nan, also stationary file will be written.
    
    """          
    
    f.writelines('SWAN ' + str(self.version) + '\n')
    it = 0
    for it,t in enumerate(self.text):
        f.writelines('$' + t + '\n')
        
    # TIME    
    if not timecoding is None:
        if type(self.t[0]) is type(datetime.datetime(1997,1,1)):
            f.writelines('TIME\n')
            #f.writelines(str(len(self.x)) + '\n') # can SWAN handle this?, that would be very convenient for pre-allocation
            f.writelines(str(min(1,len(self.x))) + '\n')

    # LOC
    if self.lon[0] == None or np.isnan(self.lon[0]):
        f.writelines('LOCATIONS\n')
        f.writelines(str(len(self.x)) + '\n')
        for i in range(len(self.x)):
            f.write("{:g} {:g}\n".format(self.x[i],self.y[i]))
        nxy = len(self.x)
    else:
        f.writelines('LONLAT\n')
        f.writelines(str(len(self.lon)) + '\n')
        for i in range(len(self.lon)):
            f.write("{:f} {:f}\n".format(self.lon[i],self.lat[i]))
        nxy = len(self.lon)            
        
    # FREQ
    try:
        if self.ftype.lower() == 'absolute':
            f.writelines('AFREQ\n')
        else:    
            f.writelines('RFREQ\n')
    except:
            f.writelines('AFREQ\n')
    f.writelines(str(len(self.f))+ '\n')
    for freq in self.f:
        f.writelines("{:g}\n".format(freq))

    # QUANT
    f.writelines('QUANT \n')
    f.writelines('3 \n')

    f.writelines('VaDens \n')
    f.writelines(self.energy_units + '\n')
    f.writelines('NaN' + '\n')

    if self.direction_units.lower() == 'degrees_true':
        f.writelines('CDIR\n')
    elif self.direction_units.lower() == 'degrees_north':
        f.writelines('NDIR\n')
    else:    
        f.writelines('NDIR\n')
        print('no units defined, assumed NDIR')
    f.writelines(self.direction_units + '\n')
    f.writelines('NaN' + '\n')

    f.writelines('DSPRDEGR \n')
    f.writelines(self.spreading_units + '\n')
    f.writelines('NaN' + '\n')


    for it,t in enumerate(self.t):
        if not timecoding is None:
            if type(self.t[0]) is type(datetime.datetime(1997,1,1)):
                if debug:
                    print('time ',self.t[it])        
                f.writelines("{:%Y%m%d.%H%M%S}                         date and time\n".format(self.t[it]))
        for ix in range(nxy):
            f.writelines('LOCATION {:d} \n'.format(ix+1))
            for ifr,fr in enumerate(self.f):
                f.writelines('{:f} {:} {:g} \n'.format(self.energy[it,ix,ifr],self.direction[it,ix,ifr],self.spreading[it,ix,ifr]))
                
    return
    
#@static
def to_file0D(self,f):
    """
    Save paramatetric tiemseries to SWAN TPAR file.
    
    S0 = ow.Spec0() 
    # fill S1
    f = fopen(r'c:\temp\newfile.tpar')
    S0.to_file1D(f)
    f.close()        
    
    A TPar file has 5 columns where we use CAPITAL choices:
    
        Time (ISO-notation), 
        Hs, 
        Period (average or PEAK PERIOD depending on the choice given in command BOUND SHAPE), 
        Peak Direction (NAUTICAL or Cartesian, depending on command SET), 
        Directional spread (in degrees or as POWER OF COS depending on the choice given in command BOUND SHAPE).
    
    """          

    # write actual meaning to file, instead of relying on input file
    # TPAR file can not have a header, but adding it after keyword TPAR works
    f.writelines('TPAR $yyyymmdd.HHMMSS Hs Tp pdir ms\n')
    if type(self.Hs) is type(1.) or type(self.Hs) is type(1):
        f.writelines(datetime.datetime.strftime(self.t,'%Y%m%d.%H%M%S')+' '+
        str(self.Hs)+' '+
        str(self.Tp)+' '+
        str(self.pdir)+' '+
        str(self.ms) + '\n') 
    else:
        for i in range(len(self.t)):
            #print(i)
            f.writelines(str(datetime.datetime.strftime(self.t[i],'%Y%m%d.%H%M%S')) + ' '+
                         str(self.Hs[i])+' '+
                         str(self.Tp[i])+' '+
                         str(self.pdir[i])+' '+
                         str(self.ms[i])+'\n')
            #datetime.datetime.strftime(self.t[i],'%Y%m%d.%H%M%S')+' '+
    
    
def from_file0D(f,source='from_swan'):
    """
    Load paramatetric tiemseries from SWAN TPAR file.
    
    f = fopen(r'c:\temp\myfile.spc')
    S0 = swan.from_swan(f,source='myfile.spc') # optionally remember source
    f.close()
    """ 
    
    self = ow.Spec0()

    raw = '$'
    while raw.strip()[0]=='$':
        raw = f.readline()
    
    t    = []
    Hs   = []
    Tp   = []
    pdir = []
    ms   = []
    raw = f.readline()
    while len(raw.strip())>0:
        t.append   (datetime.datetime.strptime(raw.split()[0],'%Y%m%d.%H%M%S'))
        Hs.append  (float(raw.split()[1]))
        Tp.append  (float(raw.split()[2]))
        pdir.append(float(raw.split()[3]))
        ms.append  (float(raw.split()[4]))
        raw = f.readline()
    f.close()

    self.t    = np.asarray(t)
    self.Hs   = np.asarray(Hs)
    self.Tp   = np.asarray(Tp)
    self.pdir = np.asarray(pdir)
    self.ms   = np.asarray(ms)    
    
    return self

if __name__ == '__main__':

    from swantest import TestSuite
    import unittest
    unittest.main()
    