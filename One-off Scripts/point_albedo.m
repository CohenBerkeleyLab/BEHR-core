% Make the albedo test files have the MODISAlbedo field point to the BRDF
% albedo
function point_albedo
p='/Volumes/share-sat/SAT/BEHR/AlbedoTest';
psave='/Volumes/share-sat/SAT/BEHR/AlbedoTest3';
F=dir(fullfile(p,'*.mat'));
for a=1:numel(F)
    fprintf('Loading %s\n',F(a).name);
    load(fullfile(p,F(a).name));
    for b=1:numel(Data)
        Data(b).MODISAlbedo = Data(b).MODISAlbedoPoly;
    end
    save(fullfile(psave,F(a).name),'Data');
end
end
