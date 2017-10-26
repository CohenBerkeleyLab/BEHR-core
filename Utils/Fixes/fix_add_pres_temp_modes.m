function [  ] = fix_add_pres_temp_modes( behr_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

F = dirff(fullfile(behr_dir, '*.mat'));
for a=1:numel(F)
    D = load(F(a).name);
    for b=1:numel(D.Data)
        D.Data(b).BEHRWRFTemperatureMode = 'precomputed';
        D.Data(b).BEHRWRFPressureMode = 'precomputed';
        D.Data(b).GitHead_WRFUtils_Main = 'N/A';
    end
    for b=1:numel(D.OMI)
        D.OMI(b).BEHRWRFTemperatureMode = 'precomputed';
        D.OMI(b).BEHRWRFPressureMode = 'precomputed';
        D.OMI(b).GitHead_WRFUtils_Main = 'N/A';
    end
    save(F(a).name, '-struct', 'D');
end

end

