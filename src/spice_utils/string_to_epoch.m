function et = string_to_epoch(tstr)
    et = cspice_str2et(convertStringsToChars(tstr));
end