function cfg = parseThemeColors(data, defaults, overrideKey)
arguments
    data struct
    defaults struct
    overrideKey char = 'themeColors'
end

cfg = defaults;
cfg = mergeStructs(cfg, data);

if isfield(data, overrideKey) && isstruct(data.(overrideKey))
    cfg = mergeStructs(cfg, data.(overrideKey));
end
end