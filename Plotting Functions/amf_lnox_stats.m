function [ stat_table_amf, stat_table_vcd ] = amf_lnox_stats(  )
%AMF_LNOX_STATS Calculates statistics about individual pixel changes when
%LNOx is included in the profile.
%   Detailed explanation goes here

hybrid_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/LNOx-AMFs/DC3-old-new-hybrid';
base_dir = BEHR_paths('behr_mat_dir');

F = dir(fullfile(hybrid_dir, 'OMI_BEHR*'));

hy_amfs = [];
hy_amfs_allcld = [];
base_amfs = [];
base_amfs_allcld = [];

hy_vcds = [];
hy_vcds_allcld = [];
base_vcds = [];
base_vcds_allcld = [];

for a=1:numel(F);
    fprintf('%s\n',F(a).name);
    hy = load(fullfile(hybrid_dir, F(a).name),'Data');
    base = load(fullfile(base_dir, F(a).name),'Data');
    
    for b=1:numel(hy.Data)
        data_hy = hy.Data(b);
        data_base = base.Data(b);
        [hy_amfs, base_amfs] = cat_valid_amfs(hy_amfs, base_amfs, data_hy, data_base, 0.2, 'BEHRAMFTrop');
        [hy_vcds, base_vcds] = cat_valid_amfs(hy_vcds, base_vcds, data_hy, data_base, 0.2, 'BEHRColumnAmountNO2Trop');
        [hy_amfs_allcld, base_amfs_allcld] = cat_valid_amfs(hy_amfs_allcld, base_amfs_allcld, data_hy, data_base, 1.0, 'BEHRAMFTrop');
        [hy_vcds_allcld, base_vcds_allcld] = cat_valid_amfs(hy_vcds_allcld, base_vcds_allcld, data_hy, data_base, 1.0, 'BEHRColumnAmountNO2Trop');
    end
end

nans = isnan(hy_amfs) | isnan(base_amfs);
hy_amfs(nans) = [];
base_amfs(nans) = [];

nans = isnan(hy_amfs_allcld) | isnan(base_amfs_allcld);
hy_amfs_allcld(nans) = [];
base_amfs_allcld(nans) = [];

nans = isnan(hy_vcds) | isnan(base_vcds);
hy_vcds(nans) = [];
base_vcds(nans) = [];

nans = isnan(hy_vcds_allcld) | isnan(base_vcds_allcld);
hy_vcds_allcld(nans) = [];
base_vcds_allcld(nans) = [];

stat_table_amf = do_all_stats(hy_amfs, hy_amfs_allcld, base_amfs, base_amfs_allcld);
stat_table_vcd = do_all_stats(hy_vcds, hy_vcds_allcld, base_vcds, base_vcds_allcld);

end

function [hy_amfs, base_amfs] = cat_valid_amfs(hy_amfs, base_amfs, data_hy, data_base, maxcld, field)
data_hy.Areaweight = ones(size(data_hy.Longitude));
data_hy = omi_pixel_reject(data_hy,'omi',maxcld,'XTrackFlags',[1 58]);
xx_hy = reshape(data_hy.Areaweight > 0, size(data_hy.Longitude));

data_base.Areaweight = ones(size(data_base.Longitude));
data_base = omi_pixel_reject(data_base,'omi',maxcld,'XTrackFlags',[1 58]);
xx_base = reshape(data_base.Areaweight > 0, size(data_base.Longitude));

xx = xx_hy & xx_base;

if sum(xx(:)) > 0
    hy_amfs = cat(1, hy_amfs, data_hy.(field)(xx));
    base_amfs = cat(1, base_amfs, data_base.(field)(xx));
end
end

function T = do_all_stats(hy, hy_allcld, base, base_allcld)
% Do statistics
perdiff = reldiff(hy, base)*100;
perdiff_allcld = reldiff(hy_allcld, base_allcld)*100;

pos = perdiff > 0;
pos_all = perdiff_allcld > 0;
neg = perdiff < 0;
neg_all = perdiff_allcld < 0;

A = [];
A = cat(1, A, do_stats(perdiff));
A = cat(1, A, do_stats(perdiff_allcld));
A = cat(1, A, do_stats(perdiff(pos)));
A = cat(1, A, do_stats(perdiff_allcld(pos_all)));
A = cat(1, A, do_stats(perdiff(neg)));
A = cat(1, A, do_stats(perdiff_allcld(neg_all)));

datasets = {'Perdiff cldfrac < 0.2', 'Perdiff any cldfrac','Positive perdiff cldfrac < 0.2','Positive perdiff any cldfrac','Negative perdiff cldfrac <= 0.2','Negative perdiff any cldfrac'};

figure; boxplot(padcat(2, perdiff, perdiff_allcld, perdiff(pos), perdiff_allcld(pos_all), perdiff(neg), perdiff_allcld(neg_all)));
set(gca,'xticklabel',datasets,'ygrid','on','yminorgrid','on','yminortick','on');
title('Percent difference in AMFs');

T = array2table(A, 'VariableNames', {'Mean', 'Std', 'Median', 'Max', 'Min', 'Percentile5th','Percentile25th','Percentile75th','Percentile95th'},...
    'RowNames', datasets);
end

function stats = do_stats(vec)
vec = vec(:);
stats = [nanmean(vec), nanstd(vec), nanmedian(vec), nanmax(vec), nanmin(vec),...
         quantile(vec, [0.05 0.25 0.75 0.95])];
end