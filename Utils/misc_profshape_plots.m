function [  ] = misc_profshape_plots( plottype )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

prof2013_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EMG fitting nationally/BEHR_100Emis';
prof2011scaled_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EmissionsScaling/With Scaling';
prof2011unscaled_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EmissionsScaling/Control';

switch lower(plottype)
    case 'emis_effect'
        plot_emis_effect()
    case 'poc_hist'
        proof_of_concept_hist();
    case 'poc_avg'
        proof_of_concept_avg();
    case 'windvar-sub12km'
        wind_variability_sub12km();
end

    function plot_emis_effect
        % This will plot the effect of different emissions on various levels of the
        % US June monthly profiles
        
        n_levels = 29;
        
        plotting_style = ask_multichoice('Plot as series in one plot, one plot per level, slopes vs level, or use plot slice gui?', {'Series','Plot per level','Slope vs level','Plot slice gui'}, 'list', true);
        oneplot_bool = strcmpi(plotting_style, 'series');
        metaplot_bool = any(strcmpi(plotting_style, {'slope vs level','plot slice gui'}));
        
        if ~metaplot_bool
            fprintf('What range of levels do you wish to plot?\n');
            bottom_lev = ask_number(sprintf('Enter the bottom level (1-%d)', n_levels), 'testfxn', @(x) isscalar(x) && x>=1  && x <= n_levels, 'testmsg', 'Must be a scalar between 1 and 29');
            top_lev = ask_number(sprintf('Enter the top level (%d-%d)', bottom_lev, n_levels), 'testfxn', @(x) isscalar(x) && x>=bottom_lev && x <= n_levels, 'testmsg', 'Must be a scalar greater than the bottom level and <= 29');
        else
            bottom_lev = 1;
            top_lev = n_levels;
        end
        
        % First load the data
        e100 = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/US_BEHR_FULLEMIS/monthly/WRF_BEHR_monthly_2013-06-30.nc');
        no2_100 = double(ncread(e100.Filename, 'no2'))*1e3;
        e50 = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/US_BEHR_HALFEMIS/monthly/WRF_BEHR_monthly_2013-06-30.nc');
        no2_50 = double(ncread(e50.Filename, 'no2'))*1e3;
        
        if ~strcmpi(plotting_style,'plot slice gui')
            colors = {'k','b','r',[0 0.5 0],[1 0.5 0], 'm'};
            markers = {'o', 's', 'x', '^', 'p', 'v'};
            
            lines = gobjects(2*(top_lev - bottom_lev + 1), 1);
            lstr = cell(2*(top_lev - bottom_lev + 1), 1);
            slopes = nan(n_levels, 1);
            slopes_sd = nan(n_levels, 1);
            figure;
            % Index for the line and legend strings
            i = 1;
            for a = bottom_lev:top_lev
                
                [fitx, fity, fitlstr, fitdata] = calc_fit_line(reshape(no2_100(:,:,a),[],1), reshape(no2_50(:,:,a),[],1), 'regression', 'rma', 'xcoord', 'dataxlim');
                if metaplot_bool
                    slopes(i) = fitdata.P(1);
                    slopes_sd(i) = fitdata.StdDevM;
                    i = i+1;
                else
                    % Index for the colors and markers. Let be one (same for every plot) if
                    % plotting multiple figures
                    if oneplot_bool
                        j = mod(a-1, numel(markers))+1;
                    else
                        j = 1;
                    end
                    
                    lines(i) = line(reshape(no2_100(:,:,a),[],1), reshape(no2_50(:,:,a),[],1), 'marker', markers{j}, 'color', colors{j}, 'linestyle','none');
                    lstr{i} = sprintf('Level %d', a);
                    
                    lines(i+1) = line(fitx, fity, 'color', colors{j}, 'linestyle', '--');
                    
                    eol = strfind(fitlstr, char(10)); % looks for newline
                    lstr{i+1} = fitlstr(1:eol-1);
                    
                    if ~oneplot_bool
                        % Make legend now if doing multiple plots; do legend after loop if
                        % putting all series on one plot
                        legend(lines(i:i+1), lstr(i:i+1));
                    end
                    
                    if ~oneplot_bool && a < top_lev
                        xlabel('[NO_2] (ppbv) 2013 emissions');
                        ylabel('[NO_2] (ppbv) 50% 2013 emissions');
                        set(gca,'fontsize',16);
                        figure;
                    end
                    i = i+2;
                end
            end
            
            if metaplot_bool
                line(bottom_lev:top_lev, slopes, 'marker', 's', 'linestyle', 'none', 'color', 'k');
                scatter_errorbars(bottom_lev:top_lev, slopes, slopes_sd, 'color', 'k');
                xlabel('Level index');
                ylabel('Slope [NO_2]_{50% emis} vs [NO_2]_{100% emis}')
                set(gca,'fontsize',16);
            else
                if oneplot_bool
                    legend(lines, lstr)
                end
                xlabel('[NO_2] (ppbv) 2013 emissions');
                ylabel('[NO_2] (ppbv) 50% 2013 emissions');
                set(gca,'fontsize',16);
            end
        else
            lon = ncread(e50.Filename,'XLONG');
            lat = ncread(e50.Filename,'XLAT');
            [slon, slat] = state_outlines('not','ak','hi');
            plot_slice_gui(no2_50 ./ no2_100, lon, lat, slon, slat,'Ratio of [NO_2] 50% emis / 100% emis');
            w = reldiff(no2_50, no2_100) ./ -0.5;
            plot_slice_gui(w, lon, lat, slon, slat,'Weight');
        end
    end

    function proof_of_concept_hist()
        % Load the three data sets, compare the differences between the
        % true 2013 daily profile retrieval and the 2011 scaled and
        % unscaled retrieval
        vcds_2013 = load_sat_data(prof2013_dir);
        vcds_2011scaled = load_sat_data(prof2011scaled_dir);
        vcds_2011unscaled = load_sat_data(prof2011unscaled_dir);
        
        del_scaled = vcds_2011scaled - vcds_2013;
        perdel_scaled = reldiff(vcds_2011scaled, vcds_2013)*100;
        del_unscaled = vcds_2011unscaled - vcds_2013;
        perdel_unscaled = reldiff(vcds_2011unscaled, vcds_2013)*100;
        
        figure; 
        histsemilog(del_scaled(:),50);
        title('Difference between using scaled 2011 profiles and 2013 profiles');
        
        figure; 
        histsemilog(perdel_scaled(:),50);
        title('Percent difference between using scaled 2011 profiles and 2013 profiles');
        
        figure; 
        histsemilog(del_unscaled(:),50);
        title('Difference between using unscaled 2011 profiles and 2013 profiles');
        
        figure; 
        histsemilog(perdel_unscaled(:),50);
        title('Percent difference between using unscaled 2011 profiles and 2013 profiles');
    end

    function proof_of_concept_avg()
        % Makes monthly average difference maps of the scaled and unscaled
        % 2011 profile derived retrievals
        [~, no2_true, lon, lat] = no2_column_map_2014('2013-06-01','2013-06-30',[-125 -65],[25 50],'behrdir',prof2013_dir);
        [~, no2_scaled] = no2_column_map_2014('2013-06-01','2013-06-30',[-125 -65],[25 50],'behrdir',prof2011scaled_dir);
        [~, no2_unscaled] = no2_column_map_2014('2013-06-01','2013-06-30',[-125 -65],[25 50],'behrdir',prof2011unscaled_dir);
        
        del_scaled = no2_scaled - no2_true;
        perdel_scaled = reldiff(no2_scaled, no2_true)*100;
        del_unscaled = no2_unscaled - no2_true;
        perdel_unscaled = reldiff(no2_unscaled, no2_true)*100;
        
        del_plot(del_scaled, lon, lat, 'Difference scaled 2011 - daily 2013 profiles');
        del_plot(perdel_scaled, lon, lat, 'Percent difference scaled 2011 - daily 2013 profiles');
        del_plot(del_unscaled, lon, lat, 'Difference unscaled 2011 - daily 2013 profiles');
        del_plot(perdel_unscaled, lon, lat, 'Percent difference unscaled 2011 - daily 2013 profiles');
    end

    function wind_variability_sub12km()
        % Using one of Xueling's 3 km runs, calculate how much variability
        % there is in wind within a single 12 km WRF grid cell.
        wrfpath='/Volumes/share2/USERS/LaughnerJ/WRF/3km/domain_setup';
        wrffile='070215_d02.nc';
        [U,V,COSALPHA,SINALPHA,xlon,xlat] = read_wrf_vars(wrfpath, {wrffile}, {'U','V','COSALPHA','SINALPHA','XLONG','XLAT'});
        [Ue, Ve] = wrf_winds_transform(U,V,COSALPHA,SINALPHA);
        
        winddir = atan2d(Ve, Ue);
        windvel = sqrt(Ue.^2 + Ve.^2);
        winddir = winddir(:,:,1);
        windvel = windvel(:,:,1);
        
        % There are 16 3 km grid cells per 1 12 km grid cell, we will take
        % the std. dev. of the 3 km winds in 4x4 chunks.
        blank_vec = nan(ceil(numel(winddir)/16),1);
        winddir_std = blank_vec;
        windvel_std = blank_vec;
        
        i=0;
        for a=1:size(winddir,1)
            a2 = min(a+3,size(winddir,1));
            for b=1:size(winddir,2)
                b2 = min(b+3,size(winddir,2));
                i=i+1;
                winddir_tmp = reshape(winddir(a:a2,b:b2),[],1);
                windvel_tmp = reshape(windvel(a:a2,b:b2),[],1);
                
                % Wind direction needs handled a little differently to
                % account for the wrapping around of angles
                winddir_std(i) = sqrt( nansum2(angle_diffd(winddir_tmp, angle_meand(winddir_tmp)).^2) ./ (numel(winddir_tmp)-1) );
                windvel_std(i) = nanstd(windvel_tmp);
            end
        end
        
        figure;
        boxplot(winddir_std);
        title('Wind direction')
        figure;
        boxplot(windvel_std);
        title('Wind speed');
    end
end

function vcds = load_sat_data(sat_dir)
F = dir(fullfile(sat_dir,'OMI_BEHR*.mat'));
vcds = [];
for a=1:numel(F)
    fprintf('Loading %s\n',F(a).name);
    D = load(fullfile(sat_dir, F(a).name),'Data');
    Data = D.Data;
    for b=1:numel(Data)
        Data(b).Areaweight = ones(size(Data(b).Longitude));
        Data(b) = omi_pixel_reject(Data(b),'omi',0.2,'XTrackFlags');
        Data(b).BEHRColumnAmountNO2Trop(Data(b).Areaweight==0) = nan;
        vcds = cat(1, vcds, Data(b).BEHRColumnAmountNO2Trop);
    end
end
end

function del_plot(del, lon, lat, titlestr)
figure; 
pcolor(lon, lat, del);
shading flat;
caxis([-max(abs(del(:))), max(abs(del(:)))]);
state_outlines('not','ak','hi');
colormap(blue_red_cmap);
title(titlestr);
end