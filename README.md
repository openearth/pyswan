# oceanwaves, applied to swan
Generic toolbox for spectral oceanwaves plus SWAN IO toolbox

The oceanwaves package contains general classes for 1D and 2D spectral wave data.
The classes come with dedicated plot functions (see notebook), and functions to integrate spectral 
parameters (Hm0, Tp) on-the-fly. The oceanwaves classes are meant to be used for a series of 
packages for IO from different wave data sources. Here we implemented the use of oceanwaves as a 
basis for swan IO, but we foresee additional packages for various bouy suppliers, metocean 
data suppliers (web services), OGC and SWE services, other wave models (xbeach, wam), 
and FMEA codes for analysis of marine structures.
