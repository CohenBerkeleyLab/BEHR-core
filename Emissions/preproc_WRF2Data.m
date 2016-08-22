function [  ] = preproc_WRF2Data(  )
%PREPROC_WRF2DATA Process WRF BEHR output into Data structures
%   The goal here is to get WRF data into a format that can be used to
%   calculate line densities as a sort of a priori-less dataset to compare
%   agains the various a priori results. This will use WRF files already
%   processes by the (slurm)run_wrf_output.sh file, since we need to be
%   able to compute NO2 columns.

wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/hourly/';
save_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US WRF Hourly';
allowed_hours = 19; % can be an array of UTC hours to include; or an empty array to use all. 
% Note that line density calculation uses all swaths (which is how UTC
% hours are saved) so using all hours will end up giving an average over
% several hours, which isn't ideal for comparison against satellite
% measurements.
start_date = '2013-06-01';
end_date = '2013-08-30';

F = dir(fullfile(wrf_path,'WRF_BEHR*.nc'));
sdnum = datenum(start_date);
ednum = datenum(end_date);
for a=1:numel(F)
    [s,e] = regexp(F(a).name,'\d\d\d\d-\d\d-\d\d');
    dnum = datenum(F(a).name(s:e));
    if dnum < sdnum || dnum > ednum
        continue
    end
    
    % Get the geographic and temporal information and the NO2 tropospheric
    % columns
    [xlon, xlat, z, utchr] = read_wrf_vars(wrf_path, F(a), {'XLONG','XLAT','z','utchr'});
    xlon = xlon(:,:,1); % these don't change with time.
    xlat = xlat(:,:,1);
    z = z(:,:,1,1);
    
    [xloncorn, xlatcorn] = wrf_grid_corners(xlon(:,:,1), xlat(:,:,1));
    trop_no2 = compute_wrf_trop_colums(fullfile(wrf_path,F(a).name));
    surfpres = 1013*exp(-z ./ 7400);
    
    % These are the fields expected by hdf_quadrangle_OMI_winds, which is
    % used to grid the Data structures after rotating the pixels. Clouds
    % and quality flags are set to 0 to indicate that all "pixels" are
    % valid. Cloud pressure is set to 400 hPa as a fill value; it shouldn't
    % be used since cloud fractions are 0. VZA needs to be 1 because a
    % later step filters out large pixels by VZA >= 60 deg and also ignores
    % swaths with all VZA = 0 (because at some point in BEHR's life that
    % meant that the swath hadn't really been processed)
    nanfill = nan(size(xlon));
    zerofill = zeros(size(xlon));
    onefill = ones(size(xlon));
    Data = struct('Latitude',xlat,'Longitude',xlon,'Latcorn',xlatcorn,'Loncorn',xloncorn,...
        'BEHRColumnAmountNO2Trop',nanfill,'ViewingZenithAngle',onefill,'SolarZenithAngle',nanfill,'AMFTrop',nanfill,'CloudFraction',zerofill,...
        'CloudRadianceFraction',zerofill,'CloudPressure',400*onefill,'ColumnAmountNO2Trop',nanfill,'RelativeAzimuthAngle',nanfill,'MODISAlbedo',nanfill,...
        'GLOBETerpres',surfpres,'BEHRAMFTrop',nanfill,'vcdQualityFlags',zerofill,'XTrackQualityFlags',zerofill,'UTC_hr',[]);
    
    u_utchrs = unique(utchr);
    u_utchrs = u_utchrs(ismember(u_utchrs, allowed_hours));
    Data = repmat(Data,numel(u_utchrs),1);
    for s=1:numel(Data)
        tt = utchr == u_utchrs(s);
        Data(s).BEHRColumnAmountNO2Trop = nanmean(trop_no2(:,:,tt),3);
        Data(s).UTC_hr = u_utchrs(s);
    end
    
    save_name = sprintf('WRF_EMGData_%04d%02d%02d.mat',year(dnum),month(dnum),day(dnum));
    fprintf('Saving %s\n',fullfile(save_path,save_name));
    save(fullfile(save_path,save_name),'Data');
end


end

