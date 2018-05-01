function [ testAMFs, testVisAMFs ] = test_published_scattering_weights( Data )
% TEST_PUBLISHED_SCATTERING_WEIGHTS Reproduces AMFs from the published scattering weights and profiles
%   TESTAMFS = TEST_PUBLISHED_SCATTERING_WEIGHTS( DATA ) Given a structure
%   DATA that contains the fields BEHRScatteringWeights, BEHRNO2apriori,
%   BEHRPressureLevels, and GLOBETerpres, this will calculate total
%   tropospheric AMFs from those values, which will be returned as
%   TESTAMFS. You can compare these to the BEHRAMFTrop field to verify that
%   the two are the same to within a few percent. (There is always some
%   difference because in BEHR proper we handle the clear and cloudy AMFs
%   separately, which reduces some interpolation errors around cloud
%   pressure where the combined scattering weights have a sharp
%   discontinuity.

if ~isscalar(Data) || ~isstruct(Data)
    error('input:bad_type', 'DATA must be a scalar structure');
end

testAMFs = nan(size(Data.BEHRAMFTrop));
testVisAMFs = nan(size(Data.BEHRAMFTrop));

for i=1:numel(Data.BEHRAMFTrop)
    notnans = ~isnan(Data.BEHRPressureLevels(:,i));
    if sum(notnans) > 1
        % If we don't have tropopause pressure defined, then we need to
        % revert to the old way of doing things where we just assume a
        % pressure of 200 hPa
        if isfield(Data, 'BEHRTropopausePressure')
            tropopausePres = Data.BEHRTropopausePressure(i);
        else
            tropopausePres = 200;
        end
        
        % Do the testing differently if we have separate clear and cloudy
        % scattering weights or just a single vector of scattering weights
        if isfield(Data, 'BEHRScatteringWeightsClear')
            % The published scattering weights already include the temperature
            % correction. Since the correction is calculated as 1 - 0.003*(T -
            % 220) in omiAmfAK2, passing a temperature profile of 220 at all
            % pressures will apply a temperature correction of 1 (i.e. no
            % correction).
            temperature = 220*ones(sum(notnans),1);
            % In BEHR, we have to restrict surface pressure to >= 1013 due to
            % the limitation of the scattering weight look up table.
            surfPres = min(Data.GLOBETerpres(i), 1013);
            cldPres = min(Data.CloudPressure(i), surfPres);
            
            
            
            [testAMFs(i), testVisAMFs(i)] = omiAmfAK2(surfPres, tropopausePres, cldPres, Data.CloudFraction(i), Data.CloudRadianceFraction(i),...
                Data.BEHRPressureLevels(notnans, i), Data.BEHRScatteringWeightsClear(notnans,i), Data.BEHRScatteringWeightsCloudy(notnans,i),...
                temperature, Data.BEHRNO2apriori(notnans,i));
        else
            tmp_amf = integPr2(Data.BEHRScatteringWeights(notnans,i) .* Data.BEHRNO2apriori(notnans,i), Data.BEHRPressureLevels(notnans,i), Data.GLOBETerpres(i), tropopausePres);
            tmp_model_ground_vcd = integPr2(Data.BEHRNO2apriori(notnans,i), Data.BEHRPressureLevels(notnans,i), Data.GLOBETerpres(i), tropopausePres);
            tmp_model_cloud_vcd = integPr2(Data.BEHRNO2apriori(notnans,i), Data.BEHRPressureLevels(notnans,i), Data.CloudPressure(i), tropopausePres);
            
            testAMFs(i) = tmp_amf / tmp_model_ground_vcd;
            
            tmp_model_vis_vcd = tmp_model_cloud_vcd * Data.CloudFraction(i) + tmp_model_ground_vcd * (1-Data.CloudFraction(i));
            testVisAMFs(i) = tmp_amf / tmp_model_vis_vcd;
        end
    end
end

end

