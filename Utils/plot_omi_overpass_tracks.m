function [ output_args ] = plot_omi_overpass_tracks( start_date, end_date, varargin )
%PLOT_OMI_OVERPASS_TRACKS Plots OMI overpass tracks
%   PLOT_OMI_OVERPASS_TRACKS( START_DATE, END_DATE ) plots global overpass
%   tracks for the date range given.
%
%   PLOT_OMI_OVERPASS_TRACKS( START_DATE, END_DATE, 'region', region)
%   plots overpass tracks for the given region. Currently only allowed
%   region is 'US'.

%%%%% INPUT PARSING %%%%%
E = JLLErrors;

p = inputParser;
p.addParameter('region','global');

p.parse(varargin{:});
pout = p.Results;

region = lower(pout.region);

allowed_regions = {'global','us'};
if ~ismember(region, allowed_regions)
    E.badinput('Parameter "region" must be one of %s',strjoin(allowed_regions,', '));
end

try
    sdnum=datenum(start_date);
    ednum=datenum(end_date);
catch err
    if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
        E.badinput('START_DATE and END_DATE must be valid date numbers or date strings')
    end
end

if isnan(sdnum) || isinf(sdnum) || isnan(ednum) || isinf(ednum)
    E.badinput('START_DATE and END_DATE must be valid date numbers or date strings')
end


%%%%% MAIN FUNCTION %%%%%

% First load the sat lat/lon for all the days
dnums = sdnum:ednum;
sat_lon = cell(size(dnums));
sat_lat = cell(size(dnums));

for d=1:numel(dnums)
    yr=year(dnums(d));
    mn=month(dnums(d));
    dy=day(dnums(d));
    
    fpath = fullfile(omno2_dir, sprintf('%04d',yr), sprintf('%02d',mn));
    filepat = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',yr,mn,dy);
    F = dir(fullfile(fpath,filepat));
    
    sat_lon_tmp = [];
    sat_lat_tmp = [];
    
    for f=1:numel(F)
        filetime = str2double(regexp(F(f).name,'(?<=t)\d\d\d\d(?=-)','match','once'));
        if strcmp(region,'us') && filetime < 1600
            % Skip nighttime overpasses
            continue
        end
        
        hi = h5info(fullfile(fpath, F(f).name));
        lon = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'SpacecraftLongitude')));
        lat = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'SpacecraftLatitude')));
        
        sat_lon_tmp = cat(1, sat_lon_tmp, lon(:));
        sat_lat_tmp = cat(1, sat_lat_tmp, lat(:));
    end
    
    sat_lon{d} = sat_lon_tmp;
    sat_lat{d} = sat_lat_tmp;
end

% Now make the plots.
l = gobjects(numel(dnums),1);
legstr = cell(1, numel(dnums));

% No way to get colormap without opening a figure
f=figure;
cmap=colormap('lines');
close(f);

C = load('coast');

if strcmp(region,'global')
    figure; 
    worldmap('world');
    linem(C.lat, C.long, 'color', 'k');
    for a=1:numel(sat_lon)
        l(a) = linem(sat_lat{a}, sat_lon{a}, 'color', cmap(a, :));
        legstr{a} = datestr(dnums(a), 'yyyy-mm-dd');
    end
elseif strcmp(region,'us')
    figure;
    usmercatorwide;
    for a=1:numel(sat_lon)
        l(a) = m_line(sat_lon{a}, sat_lat{a}, 'color', cmap(a,:));
        legstr{a} = datestr(dnums(a), 'yyyy-mm-dd');
    end
end

legend(l, legstr)

end

