function [  ] = plot_omi_overpass_tracks( start_date, end_date, varargin )
%PLOT_OMI_OVERPASS_TRACKS Plots OMI overpass tracks
%   PLOT_OMI_OVERPASS_TRACKS( START_DATE, END_DATE ) plots global overpass
%   tracks for the date range given.
%
%   PLOT_OMI_OVERPASS_TRACKS( ___, 'region', region)
%   plots overpass tracks for the given region. Currently only allowed
%   region is 'US'.
%
%   PLOT_OMI_OVERPASS_TRACKS( ___, 'shape', 'edges')
%   plots the approximate edges of the swath, rather than the center line.
%
%   PLOT_OMI_OVERPASS_TRACKS( ___, 'markswaths', true/false)
%   turn on or off marking swaths to differentiate each one on a given day.
%   Defaults to FALSE when plotting the center and TRUE when plotting the
%   edges.
%
%   PLOT_OMI_OVERPASS_TRACKS( ___, 'patch', true)
%   will plot transparent patches instead of lines to show the different
%   swaths. Will only do one day, and will always force shape to be 'edges'
%   if true.

%%%%% INPUT PARSING %%%%%
E = JLLErrors;

p = inputParser;
p.addParameter('region','global');
p.addParameter('shape','center');
p.addParameter('markswaths', []);
p.addParameter('patch', false);

p.parse(varargin{:});
pout = p.Results;

region = lower(pout.region);
plotshape = lower(pout.shape);
markswaths = pout.markswaths;
patchbool = pout.patch;

allowed_regions = {'global','us'};
if ~ismember(region, allowed_regions)
    E.badinput('Parameter "region" must be one of %s',strjoin(allowed_regions,', '));
end
allowed_shapes = {'center', 'edges'};
if ~ismember(plotshape, allowed_shapes)
    E.badinput('Parameter "shape" must be one of %s', strjoin(allowed_shapes, ', '));
end
if isempty(markswaths)
    if strcmp(plotshape, 'center')
        markswaths = false;
    elseif strcmp(plotshape, 'edges')
        markswaths = true;
    else
        E.notimplemented('default MARKSWATHS for shape=%s',plotshape);
    end
else
    if ~isscalar(markswaths) || ~isnumeric(markswaths) && ~islogical(markswaths)
        E.badinput('The parameter MARKSWATHS must be a scalar logical or numeric value');
    end
end
if ~isscalar(patchbool) || ~isnumeric(patchbool) && ~islogical(patchbool)
    E.badinput('The parameter PATCH must be a scalar logical or numeric value');
elseif strcmp(plotshape, 'center') && patchbool
    warning('Cannot plot patches if plotting center only; will draw patch objets around the swath edges')
    plotshape='edges';
end

try
    sdnum=datenum(start_date);
    ednum=datenum(end_date);
catch err
    if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
        E.badinput('START_DATE and END_DATE must be valid date numbers or date strings')
    else
        rethrow(err);
    end
end

if isnan(sdnum) || isinf(sdnum) || isnan(ednum) || isinf(ednum)
    E.badinput('START_DATE and END_DATE must be valid date numbers or date strings')
end


%%%%% MAIN FUNCTION %%%%%

% First load the sat lat/lon for all the days
dnums = sdnum:ednum;
swath_lon = cell(size(dnums));
swath_lat = cell(size(dnums));
swath_utc = cell(size(dnums));

for d=1:numel(dnums)
    yr=year(dnums(d));
    mn=month(dnums(d));
    dy=day(dnums(d));
    
    fpath = fullfile(omno2_dir, sprintf('%04d',yr), sprintf('%02d',mn));
    filepat = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',yr,mn,dy);
    F = dir(fullfile(fpath,filepat));
    
    swath_lon_tmp = cell(15,1); % there's usually 14--15 swaths per day, globally
    swath_lat_tmp = cell(15,1); 
    swath_utc_tmp = nan(15,2);
    
    for f=1:numel(F)
        filetime = str2double(regexp(F(f).name,'(?<=t)\d\d\d\d(?=-)','match','once'));
        if strcmp(region,'us') && filetime < 1600
            % Skip nighttime overpasses
            continue
        end
        
        hi = h5info(fullfile(fpath, F(f).name));
        if strcmp(plotshape,'center')
            lon = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'SpacecraftLongitude')));
            lat = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'SpacecraftLatitude')));
        elseif strcmp(plotshape, 'edges')
            lon = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude')));
            lat = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude')));
            lon = [reshape(lon(1:end-1,1),[],1); reshape(lon(end,1:end-1),[],1); reshape(lon(end:-1:2, end),[],1); reshape(lon(1,end:-1:1),[],1)];
            lat = [reshape(lat(1:end-1,1),[],1); reshape(lat(end,1:end-1),[],1); reshape(lat(end:-1:2, end),[],1); reshape(lat(1,end:-1:1),[],1)];
        else
            E.notimplemented('shape = %s', plotshape);
        end
        omi_time = omi_time_conv(double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Time'))));
        
        swath_lon_tmp{f} = lon(:);
        swath_lat_tmp{f} = lat(:);
        swath_utc_tmp(f,:) = [min(omi_time), max(omi_time)];
    end
    
    xx = iscellcontents(swath_lon_tmp,'isempty') & iscellcontents(swath_lat_tmp,'isempty');
    swath_lon_tmp(xx) = [];
    swath_lat_tmp(xx) = [];
    swath_utc_tmp(xx,:) = [];
    
    swath_lon{d} = swath_lon_tmp;
    swath_lat{d} = swath_lat_tmp;
    swath_utc{d} = swath_utc_tmp;
end

% Now make the plots.
l = gobjects(numel(dnums),1);
legstr = cell(1, numel(dnums));

% No way to get colormap without opening a figure
f=figure;
cmap=colormap('lines');
close(f);

if ~patchbool
    linestyles = {'-','--',':','-.'};
    
    init_map(region);
    for a=1:numel(swath_lon)
        for b = 1:numel(swath_lon{a})
            if markswaths
                c = mod(b,numel(linestyles))+1;
                ls = linestyles{c};
            else
                ls = '-';
            end
            l(a) = linem(swath_lat{a}{b}, swath_lon{a}{b}, 'color', cmap(a, :), 'linestyle', ls);
        end
        legstr{a} = datestr(dnums(a), 'yyyy-mm-dd');
    end
    
    legend(l, legstr)
else
    for a=1:numel(swath_lon)
        for b = 1:numel(swath_lon{a})
            init_map(region);
            patchm(flipud(swath_lat{a}{b}), flipud(swath_lon{a}{b}), 'FaceColor', cmap(a,:), 'FaceAlpha', 0.2);
            titlestr = sprintf('%s:\n%s - %s', datestr(dnums(a),'yyyy-mm-dd'), datestr(swath_utc{a}(b,1), 'HH:MM'), datestr(swath_utc{a}(b,2), 'HH:MM')); 
            title(titlestr);
        end
    end
end

end

function init_map(region)
if strcmp(region,'global')
    figure; 
    worldmap('world');
    C = load('coast');
    linem(C.lat, C.long, 'color', 'k');
elseif strcmp(region,'us')
    figure; 
    worldmap([15, 60],[-135 -55]);
    setm(gca,'mapprojection','mercator');
    state_outlines('k');
end
end