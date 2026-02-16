function gs_codes = ground_stations
    % From: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/ndosl_190716_v02.cmt

    % ESTRACK Deep Space Stations
    gs_codes.ESA_New_Norcia = 'NN1D';
    gs_codes.ESA_Cebreros = 'CB1D';
    gs_codes.ESA_Malargue = 'MG1D';

    % NASA DSN (70m Antennas)
    gs_codes.DSN_Goldstone = 'DS14';
    gs_codes.DSN_Canberra = 'DS43';
    gs_codes.DSN_Madrid = 'DS63';
end