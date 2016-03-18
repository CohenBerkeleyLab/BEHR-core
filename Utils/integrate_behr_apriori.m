function [ Data ] = integrate_behr_apriori( Data )
%DATA = INTEGRATE_BEHR_APRIORI( DATA ) Integrates the a priori profiles to give an a priori column density
%   Uses integPr2.m to calculate the column density represented by the a
%   priori profile.  Input DATA must be a BEHR data structure with fields
%   BEHRNO2apriori, BEHRPressureLevels, and GLOBETerpres. Will add the
%   field BEHRaprioriVCD with the result of this calculation.

E = JLLErrors;

if ~isstruct(Data) || ~isfield(Data,'BEHRNO2apriori') || ~isfield(Data,'BEHRPressureLevels') || ~isfield(Data,'GLOBETerpres')
    E.badinput('Data must be a structure with fields BEHRNO2apriori, BEHRPressureLevels, and GLOBETerpres')
end

n = numel(Data);
for a=1:n
    sz = size(Data(a).BEHRNO2apriori);
    Data(a).BEHRaprioriVCD = nan(sz(2:end));
    for p = 1:prod(sz(2:end))
        Data(a).BEHRaprioriVCD(p) = integPr2(Data(a).BEHRNO2apriori(:,p), Data(a).BEHRPressureLevels(:,p), Data(a).GLOBETerpres(p));
    end
end


end

