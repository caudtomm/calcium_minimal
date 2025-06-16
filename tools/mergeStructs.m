function out = mergeStructs(base, override)
    out = base;
    fields = fieldnames(override);
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(base, f)
            out.(f) = override.(f);
        end
    end
end
