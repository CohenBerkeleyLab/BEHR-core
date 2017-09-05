function [AMF_Data] = compare_AMF_values_sequence()
% compare_AMF_values_sequence
%
%   Allows the user to get two AMF values using different profiles,
%   geometry, cloud fraction, etc. over a series of days. SZA, VZA, PHI,
%   albedo, surface pressure, cloud pressure, cloud fraction will be the
%   same between them.  Profile 1 will be the WRF profile, profile 2 the
%   in-situ profile.
%
%   Josh Laughner <joshlaugh5@gmail.com> 12 Aug 2014

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% WARNING: this save figures to ~/Documents/MATLAB/Figures/tmp; it clears
% out that folder before running if set to true.
plot_profiles = true;

start_date = '07/01/2011';
end_date = '07/31/2011';

start_time = '12:00';
end_time = '15:00';
timezone = 'est';

no2field = 'NO2_MixingRatio_LIF';
cityname = 'CA';
shape = 'exp';

% Set to true to force cloud fractions to be 0.
nocloud = true;

% A vector of swath numbers to actually check; if you know which swath is
% the most direct overpass, it can be specified, or left empty to
% automatically go through all swaths in the data.

swaths = [3,4];

% Which site numbers to check; use this to restrict to sites that only have
% ground station data for instance.  Set as an empty matrix to use all
% sites.

sites_allowed = [];


% How to extrapolate to the top of the troposphere, acceptable choices are
% 'median', 'fit', 'wrf', 'wrf-scaled' and 'none'.  See docs for
% extrapolate_profile for options.
top_extrap = 'wrf-scaled';

% How to extrapolate to the surface for the in-situ profile.  Choices are
% 'median', 'fit', or 'ground'. See docs for extrapolate_profile for options.
bottom_extrap = 'median';

DEBUG_LEVEL = 0;

noGhost = 1; ak = 0;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% FILE PATHS %%%%%
%%%%%%%%%%%%%%%%%%%%%%

behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';
merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
profiles_dir = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles/';

amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

figure_dir = '/Users/Josh/Documents/MATLAB/Figures/tmp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CLEAR OUT TMP FIGURES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_profiles
    fig_files = dir(fullfile(figure_dir,'*.fig'));
    for a=1:numel(fig_files);
        delete(fullfile(figure_dir,fig_files(a).name));
    end
end

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%
dates = datenum(start_date):datenum(end_date);
oldprofilefile = '';

% Prep a matrix as a data table 
if isempty(swaths); nswaths = 4; else nswaths = numel(swaths); end
if isempty(sites_allowed); 
    nsites = 6; 
    sites_allowed = 0:20;
else
    nsites = numel(sites_allowed);
end
ndays = numel(dates);

header = {'Date','Profile #','WRF AMF', 'In-situ AMF', 'Residuals Magnitude', 'Signed Residuals', 'Extrapolation flag','SZA', 'VZA', 'Phi', 'MODIS Albedo', 'Surf. Pressure','Cloud Pressure', 'Cloud Frac.', 'Rad. Cld. Frac','Latitude','Longitude'};
dataTable = -9999*ones(ndays * nsites * nswaths, numel(header));
R=0;

for d=1:ndays
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    if DEBUG_LEVEL > 0; fprintf('Now on %s\n',curr_date); end
    if DEBUG_LEVEL > 1; fprintf('  Loading data...\n'); end
    
    % Load the merge and BEHR file
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('OMI_BEHR_omi*%s%s%s.mat',year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    behr_files = dir(fullfile(behr_dir,behr_filename));
    if numel(behr_files)==1
        load(fullfile(behr_dir,behr_files(1).name),'Data')
    elseif isempty(behr_files)
        if DEBUG_LEVEL > 1; fprintf('No BEHR file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of BEHR files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    % Load the profiles if we've changed months
    profile_file = sprintf('m%s_NO2_profile.mat',month);
    if ~strcmpi(profile_file, oldprofilefile) %If the profile file from the last loop is still in the same month, we don't need to reload it.
        load(fullfile(profiles_dir,profile_file));
        oldprofilefile = profile_file;
    end
    
    % Load the data from this day's merge file
    [prof_no2, utc, prof_pres, prof_lon, prof_lat] = remove_merge_fills(Merge,no2field,'alt','PRESSURE');
    prof_no2(prof_no2<0) = NaN;
    
    if DEBUG_LEVEL > 1; fprintf('  Identifying sites...'); end
    
    % Calculate what profile numbers have start times within the time
    % window chosen
    utcmin = local2utc(start_time,timezone);
    utcmax = local2utc(end_time,timezone);
    
    profnumfield = search_merge_fields(Merge,'prof');
    profnums = remove_merge_fills(Merge,profnumfield{1});
    u_profnums = unique(profnums(profnums ~= 0 & ~isnan(profnums)));
    profnums_in_time = false(size(u_profnums));
    for a=1:numel(u_profnums)
        xx = find(profnums == u_profnums(a),1,'first');
        if utc(xx) >= utcmin && utc(xx) <= utcmax;
            profnums_in_time(a) = true;
        end
    end
    final_profnums = u_profnums(profnums_in_time);
    if isempty(final_profnums); 
        if DEBUG_LEVEL > 1; fprintf(' None found.\n'); end
        continue
    else 
        if DEBUG_LEVEL > 1; fprintf('\n'); end
    end
    % Derive the site numbers from the profile number.  Baltimore profile
    % numbers are 4 digits, the first one corresponds to the site flag.  CA
    % and Houston profile numbers are six digits, the second 2 correspond
    % to the site flag.
    if final_profnums(1) < 9000 % assume MD
        final_sites = floor(final_profnums/1000);
    else % assume CA or TX
        final_sites = floor((final_profnums - 1e5)/1000);
    end
    
    
    % Now loop through all the swaths and all the profiles; find the
    % closest pixel and use that to get most of the data needed for AMF
    % calculations.  Interpolate WRF profiles to the profile's location.
    
    if isempty(swaths); swaths = 1:numel(Data); end
    for s=1:numel(swaths)
        if DEBUG_LEVEL > 0; fprintf('\tNow checking swath %d\n',swaths(s)); end
        % Check that the aircraft flight falls anywhere within the
        % satellite swath.  If not, skip the swath, if so use the results
        % of that test to import the data for the specific pixels we need
        % to compare.
        
        pix_lon = Data(swaths(s)).Longitude; pix_lat = Data(swaths(s)).Latitude;
        ss = pix_lon > min(prof_lon(:)) & pix_lon < max(prof_lon(:)) & pix_lat > min(prof_lat(:)) & pix_lat < max(prof_lat(:));
        
        if sum(ss(:)) == 0; 
            if DEBUG_LEVEL > 1; fprintf('\t  Swath %d does not overlap the aircraft flights\n',swaths(s)); end
            continue; 
        end
        
        pix_lon = pix_lon(ss); pix_lat = pix_lat(ss);
        sza = Data(swaths(s)).SolarZenithAngle(ss); vza = Data(swaths(s)).ViewingZenithAngle(ss);
        phi = Data(swaths(s)).RelativeAzimuthAngle(ss);  albedo = Data(swaths(s)).MODISAlbedo(ss);
        surfPres = Data(swaths(s)).GLOBETerpres(ss); cloudPres = Data(swaths(s)).CloudPressure(ss);
        cldFrac = Data(swaths(s)).CloudFraction(ss); cldRadFrac = Data(swaths(s)).CloudRadianceFraction(ss);
        
        if nocloud; cldFrac = zeros(size(cldFrac)); cldRadFrac = zeros(size(cldRadFrac)); end
        
        for p=1:numel(final_profnums)
            if DEBUG_LEVEL > 0; fprintf('\t\tNow on profile %d\n',final_profnums(p)); end
            if DEBUG_LEVEL > 1; fprintf('\t\t  Restricting to relevant pixels and sites...'); end
            % If the site number doesn't match any of the sites allowed,
            % skip this profile.
            if ~any(final_sites(p) == sites_allowed); 
                if DEBUG_LEVEL > 1; fprintf('  Site not allowed.\n'); end
                continue; 
            else
                if DEBUG_LEVEL > 1; fprintf('\n'); end
            end
            
            % Get the NO2, pressure, lat, and lon for this profile
            pp = profnums == final_profnums(p);
            thisprof_no2 = prof_no2(pp); thisprof_pres = prof_pres(pp);
            thisprof_lon = prof_lon(pp); thisprof_lat = prof_lat(pp);
            meanlon = nanmean(thisprof_lon); meanlat = nanmean(thisprof_lat);
            thisprof_utc = utc(pp);
            
            % Get the closest pixel for sza, vza, etc.
            delta = abs(pix_lon - meanlon) + abs(pix_lat - meanlat);
            [~,ii] = min(delta);
            sza_i = sza(ii); vza_i = vza(ii);
            phi_i = phi(ii); albedo_i = albedo(ii);
            surfPres_i = surfPres(ii); cloudPres_i = cloudPres(ii);
            cldFrac_i = cldFrac(ii); cldRadFrac_i = cldRadFrac(ii);
            
            % Get the WRF profile
            % Generate the grid matrices needed to interpolate the profile to the mean
            % lat/lon of the aircraft profile
            if DEBUG_LEVEL > 1; fprintf('\t\t  Interpolating WRF profile...\n'); end
            
            wrf_lonvec = PROFILE.Longitude(1,:);
            wrf_latvec = PROFILE.Latitude(:,1);
            wrf_pressures = fliplr(PROFILE.Pressure);
            wrf_no2 = flipdim(PROFILE.NO2_profile,1);
            [Y,Z,X] = meshgrid(wrf_latvec,wrf_pressures,wrf_lonvec);
            
            % Interpolate the WRF profiles to the location of the aircraft profile. The
            % order of the arguments may look funny, but it seems to work to match the
            % ordering of the NO2 profile matrix.
            wrf_profile = interp3(Y,Z,X,wrf_no2,nanmean(meanlat),wrf_pressures,nanmean(meanlon));
            
            % We want pressures to be monotonically increasing
            wrf_profile = flip(wrf_profile); wrf_pressures = flip(wrf_pressures);
            
            
            % Bin and extrapolate the in-situ profile
            if DEBUG_LEVEL > 1; fprintf('\t\t  Extrapolating in-situ profile...\n'); end
            [insitu_profile, insitu_pressures, flag] = extrapolate_profile(thisprof_no2, thisprof_pres, 'surfacePressure', surfPres_i,'top', top_extrap, 'bottom', bottom_extrap, 'city', cityname,...
                'utc',nanmean(thisprof_utc),'sitenum',final_sites(p), 'month',month,'lat',thisprof_lat,'lon',thisprof_lon,'date',curr_date,'shape',shape);
            
            
            % INPUT VALIDATION %
            
            % The profiles and pressure must be columns for the AMF calculation to
            % function correctly
            if ~iscolumn(wrf_pressures); wrf_pressures = wrf_pressures'; end
            if ~iscolumn(insitu_pressures); insitu_pressures = insitu_pressures'; end
            if ~iscolumn(wrf_profile); wrf_profile = wrf_profile'; end
            if ~iscolumn(insitu_profile); insitu_profile = insitu_profile'; end
            
            
            if surfPres_i > 1013;
                surfPres_p = 1013;
                flag = flag+4;
            else
                surfPres_p = surfPres_i; 
            end
            
            % WRF AMF %
            if DEBUG_LEVEL > 1; fprintf('\t\t  Calculating AMFs...\n'); end
            [temperature, ~] = rNmcTmp2(fileTmp, wrf_pressures, meanlon, meanlat, str2double(month));
            dAmfClr1 = rDamf2(fileDamf, wrf_pressures, sza_i, vza_i, phi_i, albedo_i, surfPres_p);
            cloudalbedo=0.8;
            dAmfCld1 = rDamf2(fileDamf, wrf_pressures, sza_i, vza_i, phi_i, cloudalbedo, cloudPres_i);
            
            wrf_amf = omiAmfAK2(surfPres_p, cloudPres_i, cldFrac_i, cldRadFrac_i, wrf_pressures, dAmfClr1, dAmfCld1, temperature, wrf_profile, wrf_profile, noGhost, ak);
            
            % In-Situ AMF %
            
            [temperature, ~] = rNmcTmp2(fileTmp, insitu_pressures, meanlon, meanlat, str2double(month));
            dAmfClr2 = rDamf2(fileDamf, insitu_pressures, sza_i, vza_i, phi_i, albedo_i, surfPres_p);
            cloudalbedo=0.8;
            dAmfCld2 = rDamf2(fileDamf, insitu_pressures, sza_i, vza_i, phi_i, cloudalbedo, cloudPres_i);
            
            insitu_amf = omiAmfAK2(surfPres_p, cloudPres_i, cldFrac_i, cldRadFrac_i, insitu_pressures, dAmfClr2, dAmfCld2, temperature, insitu_profile, insitu_profile, noGhost, ak);
            
            % Calculate an "R-squared" value representing the degree to
            % which the columns agree
            
            % First, interpolate the wrf profile to the pressures of the
            % in-situ profile
            wrf_insitu_pres = interp1(wrf_pressures, wrf_profile, insitu_pressures,'linear','extrap');
            
            % Normalize both profiles so that we're comparing shape factors
            wrf_profile_norm = wrf_insitu_pres ./ nansum(wrf_insitu_pres);
            insitu_profile_norm = insitu_profile ./ nansum(insitu_profile);
            
            % Now calculate the residuals as the square of the difference
            % between the profiles
            resid = nansum((wrf_profile_norm - insitu_profile_norm));
            resid2 = nansum(abs(wrf_profile_norm - insitu_profile_norm));
            
            
            %Save everything to the giant data table
            %'Date','Profile #','WRF AMF', 'In-situ AMF', 'Residuals Magnitude', 'Signed Residuals', 'SZA', 'VZA', 'Phi', 'MODIS Albedo', 'Surf. Pressure','Cloud Pressure', 'Cloud Frac.', 'Rad. Cld. Frac','Latitude','Longitude'
            R=R+1;
            row = [dates(d), final_profnums(p), wrf_amf, insitu_amf, resid2, resid, flag, sza_i, vza_i, phi_i, albedo_i, surfPres_i, cloudPres_i, cldFrac_i, cldRadFrac_i, meanlat, meanlon];
            dataTable(R,:) = row;
            
            if plot_profiles
                [base_insitu, base_pres] = bin_omisp_pressure(thisprof_pres, thisprof_no2);
                f1 = figure;
                l1 = line(wrf_profile_norm, insitu_pressures, 'color', 'b', 'linestyle','--','linewidth',2);
                l3 = line(insitu_profile_norm, insitu_pressures, 'color','r', 'linestyle','-.','linewidth',2);
                oldlim = get(gca,'xlim');
                delete(l3);
                l2 = line(base_insitu / nansum(insitu_profile), base_pres, 'color',[0 0.7 0], 'linestyle','-','linewidth',5);
                l3 = line(insitu_profile_norm, insitu_pressures, 'color','r', 'linestyle','-.','linewidth',2);
                
                legend([l1;l2;l3],{'WRF','In-situ','Extrapolated in-situ'},'fontsize',12);
                figname = sprintf('%s_%s_Prof_%d', cityname, curr_date, final_profnums(p));
                set(gca,'ydir','reverse');
                
                newlim = [oldlim(1) - abs(oldlim(2)-oldlim(1))/5, oldlim(2)];
                set(gca,'xlim',newlim);
                title(regexprep(figname,'_',' '));
                savefig(f1,fullfile(figure_dir,figname));
                close(f1);
            end
        end
    end
end

% Clean up the data table
dataTable = dataTable(1:R,:);

AMF_Table = array2table(dataTable,'variablenames',regexprep(header,'\W','_'));
AMF_Data.table = AMF_Table;
AMF_Data.top_extrap = top_extrap;
AMF_Data.bottom_extrap = bottom_extrap;
AMF_Data.shape = shape;
AMF_Data.nocloud = nocloud;
AMF_Data.no2field = no2field;
AMF_Data.swaths = swaths;
AMF_Data.sites = sites_allowed;

end
