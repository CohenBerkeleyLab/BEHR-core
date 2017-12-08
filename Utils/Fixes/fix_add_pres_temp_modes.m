function [  ] = fix_add_pres_temp_modes( behr_dir )
%FIX_ADD_PRES_TEMP_MODES Fix to add 3 attributes to BEHR files
%   FIX_ADD_PRES_TEMP_MODES( BEHR_DIR ) Loops over all .mat files in
%   BEHR_DIR and adds the fields BEHRWRFTemperatureMode,
%   BEHRWRFPressureMode, and GitHead_WRFUtils_Main to the Data and OMI
%   structures. The two "Mode" fields are given the value "precomputed" and
%   the GitHead field "N/A". These attributes needed to be added after
%   v3.0A had been run.

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

