function [  ] = check_city_winds( winds, wrf_path )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
req_fields = {'Name', 'Longitude', 'Latitude', 'dnums', 'utchr', 'windvel', 'winddir'};
if ~isstruct(winds) || any(~isfield(winds,req_fields))
    E.badinput('WINDS must be a structure output by calc_all_city_winds, and contain the fields\n%s',strjoin(req_fields, ', '));
end
if ~ischar(wrf_path) || ~exist(wrf_path, 'dir')
    E.badinput('WRF_PATH must be a valid directory')
end
F = dir(fullfile(wrf_path, 'WRF_BEHR_*.nc'));
if isempty(F)
    E.badinput('No WRF_BEHR files found in %s', wrf_path)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure out the dates of the WRF files. Compare with the dates in winds to
% find common dates.

file_dnums = nan(size(F));
for a=1:numel(F)
    date_str = regexp(F(a).name, '\d\d\d\d-\d\d-\d\d', 'match', 'once');
    file_dnums(a) = datenum(date_str, 'yyyy-mm-dd');
end

xx = ismember(file_dnums, winds(1).dnums);
common_dates = cellstr(datestr(file_dnums(xx)));

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

city_names = {winds.Name};
city = ask_multichoice('Which city?', city_names, 'list', true);

date_in = ask_multichoice('Which date?', common_dates, 'list', true);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the city data
city_ind = strcmpi(city, city_names);
city_date_ind = winds(city_ind).dnums == datenum(date_in);
city_lon = winds(city_ind).Longitude;
city_lat = winds(city_ind).Latitude;
city_utchr = winds(city_ind).utchr;
city_U = cosd(winds(city_ind).winddir(city_date_ind,:)) .* winds(city_ind).windvel(city_date_ind, :);
city_V = sind(winds(city_ind).winddir(city_date_ind,:)) .* winds(city_ind).windvel(city_date_ind, :);

% Get the WRF wind field
wrf_ind = file_dnums == datenum(date_in);
[xlon, xlat, U, V, COSALPHA, SINALPHA, zlev] = read_wrf_vars(wrf_path, F(wrf_ind), {'XLONG', 'XLAT', 'U', 'V', 'COSALPHA', 'SINALPHA', 'zlev'});
[U, V] = wrf_winds_transform(U,V,COSALPHA(:,:,1),SINALPHA(:,:,1));
[xx,yy] = find_square_around(xlon(:,:,1), xlat(:,:,1), city_lon, city_lat, 10);
[xxsm, yysm] = find_square_around(xlon(:,:,1), xlat(:,:,1), city_lon, city_lat, 1);

zlev = cumsum(zlev,3);
too_high = zlev > 500;
Unan = U;
Unan(too_high) = nan;
Vnan = V;
Vnan(too_high) = nan;

for h=1:numel(city_utchr)
    % Assume the first 5 layers of U and V are used for simple check, do 
    Usub = nanmean(U(xx,yy,1:5,h),3);
    Vsub = nanmean(V(xx,yy,1:5,h),3);
    
    figure;
    quiver(xlon(xx,yy,1), xlat(xx,yy,1), Usub, Vsub);
    hold on
    quiver(city_lon, city_lat, city_U(h), city_V(h), 'color', 'r')
    title(sprintf('UTCHR = %d', city_utchr(h)));
    
    Usub = nanmean(reshape(U(xxsm, yysm, 1:5, h),1,[]));
    Vsub = nanmean(reshape(V(xxsm, yysm, 1:5, h),1,[]));
    Usubnan = nanmean(reshape(Unan(xxsm, yysm, :, h),1,[]));
    Vsubnan = nanmean(reshape(Vnan(xxsm, yysm, :, h),1,[]));
    
    mag = reshape(sqrt(Usub.^2 + Vsub.^2),1,[]);
    nanmag = reshape(sqrt(Usubnan.^2 + Vsubnan.^2),1,[]);
    nanmag(isnan(nanmag)) = [];
    
    figure;
    line(1:3,[mag nanmag winds(city_ind).windvel(city_date_ind,h)],'marker', 'o', 'linestyle', 'none', 'color', 'k');
    set(gca,'xtick',1:3);
    set(gca,'xticklabel',{'Simple avg.','NaNed avg','Actual val'});
%     line(2:numel(mag)+1, mag, 'marker', 'o', 'linestyle', 'none', 'color', 'b')
%     line([2,numel(mag)+1], [nanmean(mag), nanmean(mag)]);
%     Usubitem = reshape(U(xxsm, yysm, 1:5, h), numel(xxsm)*numel(yysm), [])';
%     Vsubitem = reshape(V(xxsm, yysm, 1:5, h), numel(xxsm)*numel(yysm), [])';
%     X = repmat(1.25:1:9.25,5,1);
%     mag_subitem = sqrt(Usubitem.^2 + Vsubitem.^2);
%     line(X, mag_subitem, 'linestyle', 'none', 'marker', 'o', 'color', 'c')
%     
%     line(2:numel(nanmag)+1, nanmag, 'marker', 'x', 'linestyle', 'none', 'color', [0 0.5 0])
%     line([2,numel(nanmag)+1], [nanmean(nanmag), nanmean(nanmag)], 'color', [0 0.5 0]);
%     Usubitem = reshape(Unan(xxsm, yysm, 1:5, h), numel(xxsm)*numel(yysm), [])';
%     Vsubitem = reshape(Vnan(xxsm, yysm, 1:5, h), numel(xxsm)*numel(yysm), [])';
%     X = repmat(1.5:1:9.5,5,1);
%     mag_subitem = sqrt(Usubitem.^2 + Vsubitem.^2);
%     line(X, mag_subitem, 'linestyle', 'none', 'marker', 'x', 'color', 'g')
%     
%     line(1, winds(city_ind).windvel(city_date_ind,h), 'color', 'r', 'marker', 's', 'linestyle', 'none');
%     ylabel('Wind speed (m/s)')
%     title(sprintf('UTCHR = %d', city_utchr(h)));
end

tilefigs

end

