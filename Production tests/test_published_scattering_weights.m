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

testAMFs = nan(size(Data.BEHRAMFTrop));
testVisAMFs = nan(size(Data.BEHRAMFTrop));

for a=1:numel(testAMFs)
%notnans = ~isnan(Data.BEHRScatteringWeights(:,a));
notnans = true(size(Data.BEHRScatteringWeights(:,a)));
if sum(notnans) > 0
tmp_amf = integPr2(Data.BEHRScatteringWeights(notnans,a) .* Data.BEHRNO2apriori(notnans,a), Data.BEHRPressureLevels(notnans,a), Data.GLOBETerpres(a));
tmp_model_ground_vcd = integPr2(Data.BEHRNO2apriori(notnans,a), Data.BEHRPressureLevels(notnans,a), Data.GLOBETerpres(a));
tmp_model_cloud_vcd = integPr2(Data.BEHRNO2apriori(notnans,a), Data.BEHRPressureLevels(notnans,a), Data.CloudPressure(a));

testAMFs(a) = tmp_amf / tmp_model_ground_vcd;

tmp_model_vis_vcd = tmp_model_cloud_vcd * Data.CloudRadianceFraction(a) + tmp_model_ground_vcd * (1-Data.CloudRadianceFraction(a));
testVisAMFs(a) = tmp_amf / tmp_model_vis_vcd;
end
end

end

