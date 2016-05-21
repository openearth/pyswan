# -*- coding: utf-8 -*-
"""
Module for handling SWAN input/output, using oceanwwaves package:
>> import oceanwaves as ow
  
"""

"""
$HeadURL:  $
$Id:  $
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
import oceanwaves as ow
debug = False
   
#@static
def from_swan_header(f, source='from_swan'):
    """
    """
    
    self = {}
    
    self["source"] = source
    self["text"] = []
    self["version"] = f.readline().split()[1].strip()
    self["text"].append(f.readline().split('$')[1])
    self["text"].append(f.readline().split('$')[1])
    
    # when

    t = f.readline()
    if t.split()[0].strip() == 'TIME':
        timecoding = int(f.readline().split()[0].strip())
    else:
        timecoding = None

    # where

    xy  = f.readline().split()[0].strip()
    nxy = int(f.readline().split()[0])
    self["x"]   = []
    self["y"]   = []
    self["lon"] = []
    self["lat"] = []
    
    if debug:
        print('POINTS')
    for i in range(nxy):
        raw = f.readline()
        if debug:
            print(i,': ',raw.rstrip())
        if xy == 'LONLAT':
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
    
    if self["ftype"] == 'CDIR':
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
    
    if direction_type == 'CDIR':
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
    En4 = []
    if debug:
        print('TIMES')
    while not raw == '':
        if not timecoding == None:
            raw = f.readline()
            if raw == '':
                break
            tstr = raw.split()[0].strip()
            if timecoding==1:
                self.t.append(datetime.datetime.strptime(tstr,'%Y%m%d.%H%M%S'))
            else:
                self.t.append(tstr)
                print('timecoding unknown: ',timecoding)
            if debug:
                print(len(self.t),': ',self.t[-1])
                
            #print('UNDER CONSTRUCTION')
            #return self
        
        En3 = []
        for i in range(len(self.x)):
            loc = f.readline().split()[0]
            if loc =='FACTOR':
                factor = float(f.readline().split()[0].strip())
                if debug:
                    print('   POINT ',i,' FACTOR:',factor)
            else:
                if debug:
                    print('   POINT ',i,' FACTOR:',loc)
                
            En2 = []
            for j in range(len(self.f)):
                En1 = []
                raw = f.readline()
                if raw == '':
                    break
                else:
                    for k in range(len(self.direction)):
                        En1.append(float(raw.split()[k])*factor)
                En2.append(En1)
            En3.append(En2)
        En4.append(En3)
    En4 = np.asarray(En4)
    self.t      = np.asarray(self.t)
    self.energy = np.ma.masked_array(En4,En4==quantity_exception_values[0])
    
    return self 
    
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
    if quantity_names[1] == 'CDIR':
        self.direction_units = 'degrees_true'  # CF convention
    else:    
        self.direction_units = 'degrees_north' # CF convention
    if quantity_names[2] == 'DSPRD' or quantity_names[2] == 'DEGR':
        self.spreading_units = 'degrees' # DSPRD or DEGR
    else:    
        self.spreading_units = '1'       # DSPRP or POWER         
        
    # data

    raw = 'ini'
    En3,Th3,Sp3 = [],[],[]
    if debug:
        print('TIMES')
    while not raw == '':
        if not timecoding == None:
            raw = f.readline()
            if raw == '':
                break
            tstr = raw.split()[0].strip()
            if timecoding==1:
                self.t.append(datetime.datetime.strptime(tstr,'%Y%m%d.%H%M%S'))
            else:
                self.t.append(tstr)
                print('timecoding unknown: ',timecoding)
            if debug:
                print(len(self.t),': ',self.t[-1])
        En2,Th2,Sp2 = [],[],[]
        for i in range(len(self.x)):
            loc = f.readline()
#  TO DO: can be NODATA
            En1,Th1,Sp1 = [],[],[]
            for j in range(len(self.f)):
                raw = f.readline()
                if raw == '':
                    break
                else:
                    En1.append(float(raw.split()[0]))
                    Th1.append(float(raw.split()[1]))
                    Sp1.append(float(raw.split()[2]))
            En2.append(En1)
            Th2.append(Th1)
            Sp2.append(Sp1)
        En3.append(En2)
        Th3.append(Th2)
        Sp3.append(Sp2)
    En3 = np.asarray(En3)
    Th3 = np.asarray(Th3)
    Sp3 = np.asarray(Sp3)
    
    self.t         = np.asarray(self.t)
    self.energy    = np.ma.masked_array(En3,En3==quantity_exception_values[0])
    self.direction = np.ma.masked_array(Th3,Th3==quantity_exception_values[1])
    self.spreading = np.ma.masked_array(Sp3,Sp3==quantity_exception_values[2])
    
    return self
    
#@static
def to_file1D(self,f):
    """
    Save 1D spectral series file to SWAN spectral file/buffer.
    
    S1 = ow.Spec1() 
    # fill S1
    f = fopen(r'c:\temp\newfile.s1d')
    S1.to_file1D(f)
    f.close()        
    
    """          
    
    f.writelines('SWAN ' + str(self.version) + '\n')
    for t in self.text:
        f.writelines('$' + t)
        
    # TIME    
    if len(self.t)>0:
        f.writelines('TIME\n')
    f.writelines(str(len(self.x)) + '\n')

    # LOC
    print(self.lon[0])
    print(self.x[0])
    if self.x[0] == None:
        f.writelines('LONLAT\n')
        f.writelines(str(len(self.lon)) + '\n')
        for i in range(len(self.lon)):
            f.write("{:f} {:f}\n".format(self.lon[i],self.lat[i]))
    elif self.lon[0] == None:
        f.writelines('LOCATIONS\n')
        f.writelines(str(len(self.x)) + '\n')
        for i in range(len(self.x)):
            f.write("{:f} {:f}\n".format(self.x[i],self.y[i]))                 
        
    # FREQ
    if self.ftype == 'absolute':
        f.writelines('AFREQ\n')
    else:    
        f.writelines('FREQ\n')
    f.writelines(str(len(self.f))+ '\n')
    for freq in self.f:
        f.writelines("{:f}\n".format(freq))

    # QUANT
    f.writelines('QUANT \n')
    f.writelines('3 \n')

    f.writelines('VaDens \n')
    f.writelines(self.energy_units + '\n')
    f.writelines('NaN' + '\n')

    if self.direction_units == 'degrees_true':
        f.writelines('CDIR\n')
    else:    
        f.writelines('NDIR\n')
    f.writelines(self.direction_units + '\n')
    f.writelines('NaN' + '\n')

    f.writelines('DSPRDEGR \n')
    f.writelines(self.spreading_units + '\n')
    f.writelines('NaN' + '\n')


    for it,t in enumerate(self.t):
        f.writelines("{:%Y%m%d.%H%M%S}\n".format(self.t[it]))
        for ix,x in enumerate(self.x):
            f.writelines('LOCATION {:d} \n'.format(ix+1))
            for ifr,fr in enumerate(self.f):
                f.writelines('{:f} {:} {:f} \n'.format(self.energy[it,ix,ifr],self.direction[it,ix,ifr],self.spreading[it,ix,ifr]))
                
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

    import os.path
    
    OK = []

## TEST: SWAN 1D test files: different xy and directional cooordinate systems
   
    file  = r'./testdata/llcdirafreq1.spc'
    f  = open(file)
    T = from_file1D(f,source=os.path.basename(file))
    f.close()
    print(T,T.Hm0())
    T.plot(file + '.png')
    OK.append(np.abs(T.Hm0()[0,0]-1))

    file2 = r'./testdata/llcdirafreq1copy.spc'
    f2 = open(file2,'w')
    to_file1D(T,f2)
    f2.close()
    
    file = r'./testdata/xyndirrfreq1.spc'
    f = open(file)
    T = from_file1D(f,source=os.path.basename(file))
    f.close()
    print(T,T.Hm0())    
    T.plot(file + '.png')
    OK.append(np.abs(T.Hm0()[0,0]-1))
    
    file2 = r'./testdata/xyndirrfreq1copy.spc'
    f2 = open(file2,'w')
    to_file1D(T,f2)
    f2.close()    
    
## TEST: SWAN 2D test files

    file = r'./testdata/xyndirrfreq2.spc'
    f = open(file)    
    T = from_file2D(f,source=os.path.basename(file))
    f.close()
    print(T,T.energy.shape, T.Hm0())
    T.plot(file + '.png')
    OK.append(np.abs(T.Hm0()[0,0]-1))
    
## Wrap up units tests    
    
    if np.all(np.asarray(OK) < 1e-3):
        print('OK')
    else:
        raise ValueError('Not all unit tests OK')
    