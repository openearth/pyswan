import oceanwaves as ow, swan
import numpy as np, matplotlib.pyplot as plt
import datetime, os
import unittest

class TestSuite(unittest.TestCase):
    
    def test_swan1D(self):
        """Test that Hm0=1 after cycle of: generate > write > read,
        for STATionary and NOTstationary."""

        t     = [[datetime.datetime(2016,1,1),datetime.datetime(2016,1,2),datetime.datetime(2016,1,3)],
                 [np.nan],
                 [datetime.datetime(2016,1,1),datetime.datetime(2016,1,2),datetime.datetime(2016,1,3)],
                 [np.nan]]
        file  = [r'./testdata/1D_lltime.spc',
                 r'./testdata/1D_llstat.spc',
                 r'./testdata/1D_xytime.spc',
                 r'./testdata/1D_xystat.spc']
        x     = [[np.nan],[np.nan],[0,100] ,[0,100] ]
        y     = [[np.nan],[np.nan],[0,0]   ,[0,0]   ]
        lon   = [[0,100] ,[0,100] ,[np.nan],[np.nan]]
        lat   = [[0,0]   ,[0,0]   ,[np.nan],[np.nan]]
        
        f    = np.linspace(0.0250,1,40)

        for i in range(len(t)):
        
            Sp = ow.Spec1(f=f,t=t[i],y=y[i],x=x[i],lat=lat[i],lon=lon[i])
            Sp.from_jonswap(1,5,-90,10)
            #Sp.direction = Sp1.energy*0-90
            #Sp.spreading = Sp1.energy*0+10
            with open(file[i],'w') as fid:
                swan.to_file1D(Sp,fid)
            with open(file[i],'r') as fid:
                T = swan.from_file1D(fid)     
            T.plot(os.path.splitext(file[i])[0] + '.png')                
            print(file[i],T.Hm0()[0,0])
            self.assertTrue(np.abs(T.Hm0()[0,0]-1) < 1e-3)            
            # feed files to SWAN and check Hm0=1
            
    def test_swan2D(self):
        """Test that Hm0=1 after cycle of: generate > write > read,
        for STATionary and NOTstationary."""
    
        t     = [[datetime.datetime(2016,1,1),datetime.datetime(2016,1,2),datetime.datetime(2016,1,3)],
                 [np.nan],
                 [datetime.datetime(2016,1,1),datetime.datetime(2016,1,2),datetime.datetime(2016,1,3)],
                 [np.nan]]
        file  = [r'./testdata/2D_lltime.spc',
                 r'./testdata/2D_llstat.spc',
                 r'./testdata/2D_xytime.spc',
                 r'./testdata/2D_xystat.spc']
        x     = [[np.nan],[np.nan],[0,100] ,[0,100] ]
        y     = [[np.nan],[np.nan],[0,0]   ,[0,0]   ]
        lon   = [[0,100] ,[0,100] ,[np.nan],[np.nan]]
        lat   = [[0,0]   ,[0,0]   ,[np.nan],[np.nan]]
        
        dirs = list(np.arange(-12,13)*15)
        f    = np.linspace(0.0250,1,40)
    
        for i in range(len(t)):
        
            Sp = ow.Spec2(f=f,direction=dirs,t=t[i],y=y[i],x=x[i],lat=lat[i],lon=lon[i])
            Sp.from_jonswap(1,5,-90,10)
            with open(file[i],'w') as fid:
                swan.to_file2D(Sp,fid)
            with open(file[i],'r') as fid:
                T = swan.from_file2D(fid)     
            T.plot(os.path.splitext(file[i])[0] + '.png')                
            print(file[i],T.Hm0()[0,0])
            self.assertTrue(np.abs(T.Hm0()[0,0]-1) < 1e-3)            
            # feed files to SWAN and check Hm0=1          
            
if __name__ == '__main__':

    unittest.main()
    