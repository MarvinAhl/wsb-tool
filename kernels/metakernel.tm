KPL/MK

Description:
    naif0012.tls for leap seconds
    de432s.bsp for earth, sun, and moon positions
    gm_de440.tpc for the standard gravitational parameters GM
    earth_1962_250826_2125_combined.bpc for earth orientation data
    ndosl_190716_v02.bsp for the coordinates of a few ground stations
    ndosl_190716_v02.tf for the frames of a few ground stations

\begindata

    PATH_VALUES = ( '../kernels' )
    PATH_SYMBOLS = ( 'KERNELS' )
    KERNELS_TO_LOAD = (
                       '$KERNELS/naif0012.tls',
                       '$KERNELS/de432s.bsp',
                       '$KERNELS/gm_de440.tpc',
                       '$KERNELS/earth_1962_250826_2125_combined.bpc',
                       '$KERNELS/ndosl_190716_v02.bsp',
                       '$KERNELS/ndosl_190716_v02.tf'
                      )

\begintext

KERNELS: 
(downloaded from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/)