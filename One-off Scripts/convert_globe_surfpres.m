function convert_globe_surfpres(input_dir, output_dir)
%CONVERT_GLOBE_SURFPRES Change GLOBETerpres into GLOBETerrainHeight


F = dir(fullfile(input_dir, 'OMI_SP*.mat'));
parfor i_file = 1:numel(F)
    if exist(fullfile(output_dir, F(i_file).name), 'file')
        fprintf('%s exists already\n', fullfile(output_dir, F(i_file).name));
        continue
    end

    fprintf('Loading %s\n', F(i_file).name);
    D = load(fullfile(input_dir, F(i_file).name));
    Data = D.Data;
    % store the current fieldnames so that we can put GLOBETerrainHeight in
    % the same order as GLOBETerpres
    data_fields = fieldnames(Data);
    data_fields{strcmp(data_fields, 'GLOBETerpres')} = 'GLOBETerrainHeight';
    
    if ~isempty(Data)
        for i_orbit = 1:numel(Data)
            terpres = Data(i_orbit).GLOBETerpres;
            terheight = -7400 .* log(terpres ./ 1013.25);
            Data(i_orbit).GLOBETerrainHeight = terheight;
        end
        Data = rmfield(Data, 'GLOBETerpres');
        Data = orderfields(Data, data_fields);
    else
        % If Data is empty, then removing and reordering fields
        % won't work. This should create a 1x0 empty structure
        % with the right fields.
        Data = make_empty_struct_from_cell(data_fields);
        Data = Data(1, false);
    end
    saveme(fullfile(output_dir, F(i_file).name), Data);
end

end

function saveme(filename, Data) %#ok<*INUSD>
save(filename, 'Data');
end
