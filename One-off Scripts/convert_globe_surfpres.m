function convert_globe_surfpres(input_dir, output_dir)
%CONVERT_GLOBE_SURFPRES Change GLOBETerpres into GLOBETerrainHeight


F = dir(input_dir, 'OMI_SP*.mat');
parfor i_file = 1:numel(F)
    D = load(fullfile(input_dir, F(i_file).name));
    Data = D.Data;
    % store the current fieldnames so that we can put GLOBETerrainHeight in
    % the same order as GLOBETerpres
    data_fields = fieldnames(Data);
    data_fields{strcmp(data_fields, 'GLOBETerpres')} = 'GLOBETerrainHeight';
    
    for i_orbit = 1:numel(Data)
        terpres = Data(i_orbit).GLOBETerpres;
        terheight = -7400 .* log(terpres ./ 1013.25);
        Data(i_orbit).GLOBETerrainHeight = terheight;
    end
    Data = rmfield(Data, 'GLOBETerpres');
    Data = orderfields(Data, data_fields);
    saveme(fullfile(output_dir, F(i_file).name), Data);
end

end

function saveme(filename, Data) %#ok<*INUSD>
save(filename, 'Data');
end