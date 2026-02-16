function tstr = epoch_to_string(et)
    tstr = cspice_timout(et, 'YYYY-MM-DD HR:MN:SC.### ::TDB');
end