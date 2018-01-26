function [ Delta, DeltaGrid ] = behr_uncertainty_estimation( Data, OMI, parameter, percent_change )
%BEHR_UNCERTAINTY_ESTIMATION Estimate the uncertainty in BEHR NO2
%   [ DELTA, DELTAGRID ] = BEHR_UNCERTAINTY_ESTIMATION( DATA, PARAMETER, PERCENT_CHANGE )
%   This function will run the BEHR retrieval for the structure DATA but
%   with the field PARAMETER changed by PERCENT_CHANGE percent. PARAMETER
%   must match a field in DATA. This will return structures DELTA and
%   DELTAGRID which are the equivalent of Data and OMI except they will
%   have the modified parameter values and the resultant different VCDs and
%   AMFs.

E = JLLErrors;

if ~isstruct(Data)
    E.badinput('DATA must be a structure')
end
if ~ischar(parameter) || ~isfield(Data, parameter)
    % Modify later to account for the special "profile" options
    E.badinput('PARAMETER must be a field name in DATA')
end
if ~isnumeric(percent_change) || ~isscalar(percent_change)
    E.badinput('PERCENT_CHANGE must be a scalar number');
end

Delta = Data;

for a=1:numel(Delta)
    % Vary the specified parameter. Later we'll add an ability to vary the
    % NO2 profiles in a realistic way, but for now we'll stick to just
    % varying 2D parameters
    Delta(a).(parameter) = Delta(a).(parameter) * (100 + percent_change)/100;
end

% Now run BEHR but for the modified parameters
[Delta, DeltaGrid] = BEHR_main_one_day(Delta, 'profile_mode', Delta(1).BEHRProfileMode, 'lookup_profile', false);


% Calculate the percent differences in the NO2 columns and AMFs
for a=1:numel(Delta)
    Delta(a).PercentChangeNO2 = reldiff(Delta(a).BEHRColumnAmountNO2Trop, Data(a).BEHRColumnAmountNO2Trop)*100;
    Delta(a).PercentChangeNO2Vis = reldiff(Delta(a).BEHRColumnAmountNO2TropVisOnly, Data(a).BEHRColumnAmountNO2TropVisOnly)*100;
    Delta(a).PercentChangeAMF = reldiff(Delta(a).BEHRAMFTrop, Data(a).BEHRAMFTrop)*100;
    Delta(a).PercentChangeAMFVis = reldiff(Delta(a).BEHRAMFTropVisOnly, Data(a).BEHRAMFTropVisOnly)*100;
    
    DeltaGrid(a).PercentChangeNO2 = reldiff(DeltaGrid(a).BEHRColumnAmountNO2Trop, OMI(a).BEHRColumnAmountNO2Trop)*100;
    DeltaGrid(a).PercentChangeNO2Vis = reldiff(DeltaGrid(a).BEHRColumnAmountNO2TropVisOnly, OMI(a).BEHRColumnAmountNO2TropVisOnly)*100;
    DeltaGrid(a).PercentChangeAMF = reldiff(DeltaGrid(a).BEHRAMFTrop, OMI(a).BEHRAMFTrop)*100;
    DeltaGrid(a).PercentChangeAMFVis = reldiff(DeltaGrid(a).BEHRAMFTropVisOnly, OMI(a).BEHRAMFTropVisOnly)*100;
end
end

