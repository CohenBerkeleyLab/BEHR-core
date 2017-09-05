classdef misc_alb_plots
    %MISC_ALB_PLOTS Miscellaneous plots for checking the implementation of
    %the BRDF albedo in BEHR
    
    properties(Constant)
        mydir = fileparts(mfilename('fullpath'));
        black_sky_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BlackSky';
        brdf_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BRDF';
        brdf_flip_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BRDF_tomradRAA';
        workspace_dir = fullfile(behr_repo_dir, 'Workspaces', 'BRDF Albedo');
        
        scia_modis_dir = '/Users/Josh/Documents/Fortran/Sciatran/Utils/LUTGen/Tests/modis';
    end
    
    methods(Static = true)
        function make_black_sky_behr
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 but the
            % old MCD43C3 albedo.
            data_dir = misc_alb_plots.black_sky_dir;
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'fb5f363', 'fb5f363'); % commit where I'd updated to SPv3 but not changed the albedo, then cherry picked the corner fixes
            G.checkState();
            
            read_omno2_v_aug2012('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir);
            BEHR_main('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_brdf_behr
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF.
            data_dir = misc_alb_plots.brdf_dir;
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'f88f57c1', 'f88f57c1'); % commit where I'd finished merging the BRDF branch
            G.checkState();
            
            read_omno2_v_aug2012('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir);
            BEHR_main('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_soas_behr_bs
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF for the 
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'fb5f363', 'fb5f363'); % the commit where I'd fixed some of the corner calculations, but not implemented the BRDF albedo
            G.checkState();
            
            data_dir = misc_alb_plots.black_sky_dir;
            read_omno2_v_aug2012('start', '2013-05-31', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'modis_mcd43_dir', '/Volumes/share-sat/SAT/MODIS/MCD43C3');
            BEHR_main('start', '2013-05-031', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_soas_behr_brdf
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF for the 
            G = GitChecker;
            G.Strict = true;
            G.addReqCommits(behr_repo_dir, '86b60e9'); % commit where I'd made the BRDF calculation use the correct RAA definition (0 = backscattering)
            G.checkState();
            
            data_dir = misc_alb_plots.brdf_dir;
            read_omno2_v_aug2012('start', '2013-05-31', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'modis_mcd43_dir', '/Volumes/share-sat/SAT/MODIS/MCD43C1');
            BEHR_main('start', '2013-05-031', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_behr_avgs
            save_dir = misc_alb_plots.workspace_dir;
            sdate = '2005-04-01';
            edate = '2005-09-30';
            lonlim = [-125 -65];
            latlim = [25 50];
            loaddir = misc_alb_plots.black_sky_dir;
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [0 29], 'makefig', false);
            save(fullfile(save_dir, 'BlackSkyLeft.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [30 59], 'makefig', false);
            save(fullfile(save_dir, 'BlackSkyRight.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            
            loaddir = misc_alb_plots.brdf_dir;
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [0 29], 'makefig', false);
            save(fullfile(save_dir, 'BRDFLeft.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [30 59], 'makefig', false);
            save(fullfile(save_dir, 'BRDFRight.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
        end
        
        function test_raa_definition
            % Typically, the relative azimuth angle is defined as the
            % absolute value of the difference between the solar and
            % viewing azimuth angles; it is usually constrained to be
            % between 0 and 180 deg. This means that backscattering should
            % be defined as RAA = 0 and forward scattering as RAA = 180,
            % especially since Roujean et al. 1992 (which sets up the
            % kernels used in the MODIS BRDF) says on p. 20457 (last
            % paragraph first column) "...the backscattering
            % direction...phi = 0" (phi is represents the RAA).
            %
            % However, in the TOMRAD lookup table, the definition of RAA is
            % flipped so that RAA = 180 is backscattering (I believe), so
            % in BEHR our usual definition of RAA matches this. Therefore,
            % I want to test that flipping the RAA back is the correct
            % definition for the BRDF. I'm hypothesizing that if we only
            % average elements in the backscattering direction (RAA_BEHR >
            % 90) that the BRDF-derived albedo should be greater when the
            % definition passed to the BRDF kernels is s.t. RAA = 0 is
            % backscattering.
            if ask_yn('Generate the albedo averages?')
                [back_eq_0, lon, lat] = manual_average(misc_alb_plots.brdf_dir, 'MODISAlbedo', @(raa) raa > 90);
                back_eq_180 = manual_average(misc_alb_plots.brdf_flip_dir, 'MODISAlbedo', @(raa) raa > 90);
                save('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/BRDF Albedo/BRDF Implementation/raa_definition.mat', 'back_eq_0', 'back_eq_180', 'lon', 'lat');
            else
                D = load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/BRDF Albedo/BRDF Implementation/raa_definition.mat');
                back_eq_0 = D.back_eq_0;
                back_eq_180 = D.back_eq_180;
                lon = D.lon;
                lat = D.lat;
            end
            
            figure; 
            pcolor(lon,lat,back_eq_0);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 0')
            
            figure; 
            pcolor(lon,lat,back_eq_180);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 180')
            
            figure; 
            pcolor(lon,lat,back_eq_0-back_eq_180);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 0 - phi = 180')
        end
        
        function soas_scatterplots
            % Makes scatter plots of BEHR with BRDF and black sky albedo
            % vs. aircraft profiles from SOAS
            
            G = GitChecker;
            G.Strict = true;
            G.addReqCommits(no2_prof_repo, 'f52e3f12');
            G.checkState();
            
            common_opts = {'campaign', 'soas',...
                           'profnums', 'all',...
                           'profile_input', 'ranges',...
                           'merge_dir', '',...
                           'behr_prefix', 'OMI_BEHR_v2-1C_',...
                           'startdate', '',...
                           'enddate', '',...
                           'starttime', '12:00',...
                           'endtime', '15:00',...
                           'timezone', 'auto',...
                           'no2field', 'cl',...
                           'no2conv', 1e-9,...
                           'altfield', '',...
                           'radarfield', '',...
                           'presfield', '',...
                           'tempfield', '',...
                           'minheight', 0,...
                           'numBLpts', 20,...
                           'minagl', 0.5,...
                           'useground', 0,...
                           'cloudtype', 'omi',...
                           'cloudfrac', 0.2,...
                           'rowanomaly', 'XTrackFlags',...
                           'behrfield', 'BEHRColumnAmountNO2Trop',...
                           'debug', 0,...
                           'clean', 1};
            
            [prof_lon, prof_lat, ~, brdf_behr_no2, brdf_air_no2, brdf_db] = Run_Spiral_Verification('behr_dir', misc_alb_plots.brdf_dir, common_opts{:});
            [~, ~, ~, bs_behr_no2, bs_air_no2, bs_db] = Run_Spiral_Verification('behr_dir', misc_alb_plots.black_sky_dir, common_opts{:});
            
            % Scatter plot
            [brdf_fit_x, brdf_fit_y, brdf_legstr] = calc_fit_line(brdf_air_no2, brdf_behr_no2, 'regression', 'rma', 'xcoord', 0:2e15:2e16);
            [bs_fit_x, bs_fit_y, bs_legstr] = calc_fit_line(bs_air_no2, bs_behr_no2, 'regression', 'rma', 'xcoord', 0:2e15:2e16);
            l = gobjects(4,1);
            figure;
            l(1) = plot(bs_air_no2, bs_behr_no2, 'bo', 'linewidth', 2);
            hold on
            l(2) = plot(bs_fit_x, bs_fit_y, 'b--', 'linewidth', 2);
            
            l(3) = plot(brdf_air_no2, brdf_behr_no2, 'ro', 'linewidth', 2);
            l(4) = plot(brdf_fit_x, brdf_fit_y, 'r--', 'linewidth', 2);
            
            xylims;
            xlabel('Aircraft NO_2');
            ylabel('BEHR NO_2');
            set(gca,'fontsize',16);
            legend(l, {'Black sky', bs_legstr, 'BRDF', brdf_legstr});
            
            % Scatter plot with lines through origin
            [brdf_fit_x, brdf_fit_y, brdf_legstr] = calc_fit_line(brdf_air_no2, brdf_behr_no2, 'regression', 'orth-origin', 'xcoord', 0:2e15:2e16);
            [bs_fit_x, bs_fit_y, bs_legstr] = calc_fit_line(bs_air_no2, bs_behr_no2, 'regression', 'orth-origin', 'xcoord', 0:2e15:2e16);
            l = gobjects(4,1);
            figure;
            l(1) = plot(bs_air_no2, bs_behr_no2, 'bo', 'linewidth', 2);
            hold on
            l(2) = plot(bs_fit_x, bs_fit_y, 'b--', 'linewidth', 2);
            
            l(3) = plot(brdf_air_no2, brdf_behr_no2, 'ro', 'linewidth', 2);
            l(4) = plot(brdf_fit_x, brdf_fit_y, 'r--', 'linewidth', 2);
            
            xylims;
            xlabel('Aircraft NO_2');
            ylabel('BEHR NO_2');
            set(gca,'fontsize',16);
            legend(l, {'Black sky', bs_legstr, 'BRDF', brdf_legstr});
            
            % Map plots
            xl = [min(prof_lon) - 2, max(prof_lon) + 2];
            yl = [min(prof_lat) - 2, max(prof_lat) + 2];
            
            figure; scatter(prof_lon, prof_lat, [], brdf_air_no2 - brdf_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (Aircraft - BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl); 
            
            figure; scatter(prof_lon, prof_lat, [], brdf_air_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (Aircraft)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
            
            figure; scatter(prof_lon, prof_lat, [], brdf_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
            
            figure; scatter(prof_lon, prof_lat, [], bs_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
        end
        
        function varargout = ler_vs_simple_brdf(colorby, aerosol_type, DEBUG_LEVEL)
            allowed_colorbys = {'Longitude', 'Latitude', 'Site', 'SiteType', 'Month', 'SZA', 'VZA', 'RAA', 'None'};
            allowed_aer_types = {'continental', 'maritime'};
            
            if ~exist('colorby','var')
                colorby = ask_multichoice('Color the scatter plot by a variable?', allowed_colorbys, 'list', true);
            elseif ~ismember(colorby, allowed_colorbys)
                error('misc_alb_plots:bad_input', 'COLORBY must be one of %s', strjoin(allowed_colorbys, ', '));
            end
            if ~exist('aerosol_type', 'var')
                aerosol_type = ask_multichoice('Which type of aerosols for the bottom layer?', allowed_aer_types, 'list', true, 'default', 'continental');
            elseif ~ismember(aerosol_type, allowed_aer_types)
                error('misc_alb_plots:bad_input', 'AEROSOL_TYPE must be one of %s', strjoin(allowed_aer_types, ', '));
            end
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 1;
            end
            
            % Get the site output files
            aerosol_subdir = sprintf('aerosols_%s', aerosol_type);
            locs_info = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, 'loc_coeffs.nc'));
            loc_names = cellstr(ncread(locs_info.Filename, 'ShortName'));
            loc_coeffs = ncread(locs_info.Filename, 'ModisCoefficients');
            loc_months = double(ncread(locs_info.Filename, 'Months'));
                        
            for a=1:numel(loc_names)
                if DEBUG_LEVEL > 0
                    fprintf('Working on %s (%d of %d)\n', loc_names{a}, a, numel(loc_names));
                end
                loc_file = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, aerosol_subdir, sprintf('%s-LER.nc', loc_names{a})));
                loc_lers = ncread(loc_file.Filename, 'LER');
                loc_lers = permute(loc_lers, [4,3,2,1]);
                
                if a==1
                    szas = ncread(loc_file.Filename, 'sza');
                    vzas = ncread(loc_file.Filename, 'vza');
                    raas = 180 - ncread(loc_file.Filename, 'raa'); % SCIATRAN RAA definition is backwards from MODIS
                    geometries = combvec(szas', vzas', raas');
                    modis_albs = nan(numel(loc_names), numel(loc_months), size(geometries,2));
                    modis_lers = nan(size(modis_albs));
                    color_vals = nan(size(modis_albs));
                end
                
                
                for b=1:numel(loc_months)
                    this_coeffs = loc_coeffs(:,3,b,a); % always use band 3
                    modis_albs(a,b,:) = modis_brdf_alb(this_coeffs(1), this_coeffs(2), this_coeffs(3), geometries(1,:), geometries(2,:), geometries(3,:));
                    modis_lers(a,b,:) = reshape(loc_lers(:,:,:,b),[],1);
                    switch lower(colorby)
                        case 'none'
                            color_vals(a,b,:) = 0;
                        case 'site'
                            color_vals(a,b,:) = a;
                        case 'sitetype'
                            site_tmp = find(strcmp(ncreadatt(loc_file.Filename, '/', 'sitetype'), {'Cities','PowerPlants','RuralAreas'}));
                            if isempty(site_tmp)
                                error('misc_alb_plots:site_type', 'Cannot identify site type for %s', loc_names{a})
                            end
                            color_vals(a,b,:) = site_tmp;
                        case 'month'
                            color_vals(a,b,:) = loc_months(b);
                        case 'sza'
                            color_vals(a,b,:) = geometries(1,:);
                        case 'vza'
                            color_vals(a,b,:) = geometries(2,:);
                        case 'raa'
                            color_vals(a,b,:) = geometries(3,:);
                        case {'longitude', 'latitude'}
                            color_vals(a,b,:) = ncreadatt(loc_file.Filename, '/', lower(colorby));
                        otherwise
                            error('misc_alb_plots:not_implemented', 'Color by %s not implemented', colorby);
                    end
                end
            end
            
            f1 = figure;
            scatter(modis_albs(:), modis_lers(:), [], color_vals(:));
            if ~strcmpi(colorby, 'None')
                cb = colorbar;
                cb.Label.String = colorby;
            end
            xlabel('BRDF');
            ylabel('LER');
            set(gca,'fontsize',16);
            plot_fit_line(gca,[],'regression','rma');
            f1.Children(1).Location = 'northwest';
            
            perdel = reldiff(modis_lers(:), modis_albs(:), 'avg')*100;
            f2 = figure; boxplot(perdel);
            ylabel('%\Delta: LER - BRDF albedo');
            set(gca,'xtick',[],'fontsize',16, 'ygrid','on');
            if nargout > 0
                varargout = {f1, f2};
            end
        end
    end
    
end

function [avg_value, lon, lat] = manual_average(filepath, avg_field, raa_rule)
p = inputParser;
p.addParameter('avg_field', 'BEHRColumnAmountNO2Trop');
p.addParameter('raa_rule', @(raa) true(size(raa)) );

F = dir(fullfile(filepath,'OMI_BEHR_*.mat'));
do_init = true;
for a=1:numel(F)
    fprintf('Averaging %s\n', F(a).name);
    D = load(fullfile(filepath, F(a).name), 'OMI');
    if do_init
        sum_value = nan(size(D.OMI(1).Longitude));
        sum_weight = nan(size(D.OMI(1).Latitude));
        lon = D.OMI(1).Longitude;
        lat = D.OMI(1).Latitude;
        do_init = false;
    end
    for b=1:numel(D.OMI)
        omi = omi_pixel_reject(D.OMI(b),'omi',0.2,'XTrackFlags');
        
        rr = raa_rule(omi.RelativeAzimuthAngle);
        omi.Areaweight(~rr) = 0;
        sum_value = nansum2(cat(3, sum_value, omi.(avg_field) .* omi.Areaweight),3);
        sum_weight = nansum2(cat(3, sum_weight, omi.Areaweight),3);
    end
end

avg_value = sum_value ./ sum_weight;

end

