function field = unpackStruct (structure)
    fn = fieldnames(structure);
    
    for i = 1:numel(fn)
        fni = string(fn(i));
        field = structure.(fni);
        if (isstruct(field))
            unpackStruct(field);
            continue;
        end
        %assignin('base', fni, field);
    end
end