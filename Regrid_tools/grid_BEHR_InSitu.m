function [  ] = grid_BEHR_InSitu(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

DEBUG_LEVEL = 2;

start_date = '2011-07-01';
end_date = '2015-01-01';

load_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR_REPROCESSED';
load_pat = 'OMI_BEHR_*%s%s%s.mat';

save_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR_REPROCESSED';
save_prefix = 'OMI_BEHR_InSitu_gridded_';


lonmin = -125;  lonmax = -65;
latmin = 25;   latmax = 50;
resolution = 0.05; resolution2 = 0.05;

if lonmin > lonmax %Just in case I enter something backwards...
    error(E.badinput('Lonmin is greater than lonmax'))
elseif latmin > latmax
    error(E.badinput('Latmin is greater than latmax'))
end
        
dates = datenum(start_date):datenum(end_date);
for d=1:numel(dates)
    
    D = wildcard_load(load_dir, load_pat, dates(d));
    if isempty(D)
        if DEBUG_LEVEL > 1; fprintf('No file for %s\n', datestr(dates(d),'yyyy-mm-dd')); end
        continue
    end
    if DEBUG_LEVEL > 0; fprintf('Now gridding data for %s\n',datestr(dates(d),'yyyy-mm-dd')); end
    Data = D.Data;
    clear('D');
    
    if DEBUG_LEVEL > 0; disp('  Preparing OMI structure'); end
    s=numel(Data);
    
    % Prepare the OMI data structure which will receive the gridded
    % data - this will be passed to the gridding functions to keep the
    % field names in the right order.
    OMI=struct('Date','','Longitude', [], 'Latitude', [], 'Time', [], 'ViewingZenithAngle', [], 'SolarZenithAngle', [], 'ViewingAzimuthAngle', [], 'SolarAzimuthAngle', [],...
        'RelativeAzimuthAngle', [], 'AMFStrat', [], 'AMFTrop',[], 'CloudFraction', [], 'CloudRadianceFraction', [], 'CloudPressure', [], 'ColumnAmountNO2', [],...
        'SlantColumnAmountNO2', [], 'ColumnAmountNO2Trop', [], 'ColumnAmountNO2TropStd',[],'ColumnAmountNO2Strat',[],'TerrainHeight', [], 'TerrainPressure', [], 'TerrainReflectivity', [], 'vcdQualityFlags',{{}},...
        'MODISCloud', [], 'MODISAlbedo', [], 'GLOBETerpres', [], 'XTrackQualityFlags', {{}}, 'Row', [], 'Swath', [], 'TropopausePressure', [], 'BEHRColumnAmountNO2Trop',[],...
        'BEHRAMFTrop', [], 'InSituAMF', [], 'BEHR_R_ColumnAmountNO2Trop', [], 'Count', [], 'Area', [], 'Areaweight', [], 'MapData', struct);
    % Matlab treats structures as matrices, so we can duplicate our
    % structure to have the required number of entries just like a
    % matrix.
    OMI = repmat(OMI,1,s);
    hh=0;
    for s_ind=1:s;
        if Data(s_ind).ViewingZenithAngle==0;
        elseif numel(Data(s_ind).ViewingZenithAngle)==1;
            continue
        else
            if DEBUG_LEVEL > 1; fprintf('   Gridding data for swath %u\n',s_ind); end
            hh=hh+1;
            OMI(hh) = add2grid_BEHR_InSitu(Data(s_ind),OMI(hh),resolution,resolution2,[lonmin, lonmax],[latmin, latmax]); %JLL 20 Mar 2014: Superimpose data to a grid determined by lat & lon min/max and resolution above. Default resolution is 0.05 degree
        end
    end
    
    % Clean up any unused elements in OMI
    OMI(hh+1:end) = []; %#ok<NASGU>
    
    savename = sprintf('%s%04d%02d%02d.mat',save_prefix,year(dates(d)),month(dates(d)),day(dates(d)));
    if DEBUG_LEVEL > 0; disp(['   Saving data as',fullfile(save_dir,savename)]); end
    save(fullfile(save_dir,savename),'Data','OMI')
end
end

