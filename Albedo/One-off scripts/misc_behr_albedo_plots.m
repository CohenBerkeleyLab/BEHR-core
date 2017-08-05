function [  ] = misc_behr_albedo_plots( plttype )
%MISC_BEHR_ALBEDO_PLOTS Plots for the difference between black sky and BRDF

E = JLLErrors;
brdf_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF';
switch lower(plttype)
    case 'alb-quantiles'
        albedo_quantiles();
    case 'vcd-quantiles'
        vcd_quantiles()
    case 'alb-hist'
        albedo_histogram();
    case 'vcd-hist'
        vcd_histogram();
    case 'discover-validation';
        all_discover_validation();
    otherwise
        E.badinput('%s not a recognized plot type', plttype);
end



    function albedo_quantiles
        % Makes a plot of the 5th, 50th (median) and 95th percentile
        % differences in albedo for summer 2005 
        diff_type = ask_multichoice('Absolute or percent difference?',{'abs','per'});

        OMI_BRDF = cat_sat_data(brdf_dir, {'MODISAlbedoBRDF'},'prefix','OMI_BEHR','startdate','2005-06-01','enddate','2005-06-30','newdim',true,'varname','OMI','DEBUG_LEVEL',1);
        OMI_blacksky = cat_sat_data(behr_paths.behr_mat_dir, {'MODISAlbedo'},'startdate','2005-06-01','enddate','2005-06-30','newdim',true,'varname','OMI','DEBUG_LEVEL',1);
        % Load one to get the latitude/longitude
        F = dir(fullfile(brdf_dir,'OMI*.mat'));
        D = load(fullfile(brdf_dir, F(1).name),'OMI');
        lon = D.OMI(1).Longitude;
        lat = D.OMI(1).Latitude;
        clear D
        
        if strcmpi(diff_type,'abs')
            del = OMI_BRDF - OMI_blacksky;
            cblabel='BRDF - black sky, %d percentile';
        elseif strcmpi(diff_type,'per')
            del = reldiff(OMI_BRDF, OMI_blacksky)*100;
            cblabel='%%\\Delta BRDF - black sky, %d percentile';
        else
            E.notimplemented(diff_type)
        end
        q = [0.05, 0.5, 0.95];
        clim = nan(1,2);
        for a=1:numel(q)
            Z = quantile(del,q(a),3);
            fig(a)=figure; 
            pcolor(lon,lat,Z);
            shading flat
            clim(1) = min(clim(1), -max(abs(Z(:))));
            clim(2) = max(clim(2), max(abs(Z(:))));
            cb=colorbar;
            cb.Label.String = sprintf(cblabel,q(a)*100);
            set(gca,'fontsize',16);
            colormap(blue_red_cmap);
            state_outlines('k','not','ak','hi');
        end
        
        for a=1:numel(fig)
            figure(fig(a));
            caxis(clim);
        end
    end

    function albedo_histogram()
        diff_type = ask_multichoice('Absolute or percent difference?',{'abs','per'});

        OMI_BRDF = cat_sat_data(brdf_dir, {'MODISAlbedoBRDF'},'prefix','OMI_BEHR','startdate','2005-06-01','enddate','2005-06-30','varname','Data','DEBUG_LEVEL',1);
        OMI_blacksky = cat_sat_data(behr_paths.behr_mat_dir, {'MODISAlbedo'},'startdate','2005-06-01','enddate','2005-06-30','varname','Data','DEBUG_LEVEL',1);

        if strcmpi(diff_type,'abs')
            del = OMI_BRDF - OMI_blacksky;
            xstr='BRDF - black sky';
        elseif strcmpi(diff_type,'per')
            del = reldiff(OMI_BRDF, OMI_blacksky)*100;
            xstr='%\Delta BRDF - black sky';
        else
            E.notimplemented(diff_type)
        end

        figure;
        histsemilog(del(:),50);
        xlabel(xstr);
        set(gca,'fontsize',16);
    end

    function vcd_quantiles()
        diff_type = ask_multichoice('Absolute or percent difference?',{'abs','per'});
        del = load_vcd_diff('2005-06-01','2005-06-30',diff_type,'grid');
        if strcmpi(diff_type,'abs')
            cblabel='BRDF - black sky, %d percentile';
        elseif strcmpi(diff_type,'per')
            cblabel='%%\\Delta BRDF - black sky, %d percentile';
        else
            E.notimplemented(diff_type)
        end

        O = load(fullfile(brdf_dir, 'OMI_BEHR_v2-1B_20050401.mat'),'OMI');
        lon = O.OMI.Longitude;
        lat = O.OMI.Latitude;

        q = [0.05, 0.5, 0.95];
        clim = nan(1,2);
        for a=1:numel(q)
            Z = quantile(del,q(a),3);
            fig(a)=figure; 
            pcolor(lon,lat,Z);
            shading flat
            clim(1) = min(clim(1), -max(abs(Z(:))));
            clim(2) = max(clim(2), max(abs(Z(:))));
            cb=colorbar;
            cb.Label.String = sprintf(cblabel,q(a)*100);
            set(gca,'fontsize',16);
            colormap(blue_red_cmap);
            state_outlines('k','not','ak','hi');
        end
        
        for a=1:numel(fig)
            figure(fig(a));
            caxis(clim);
        end
    end

    function vcd_histogram()
        diff_type = ask_multichoice('Absolute or percent difference?',{'abs','per'});
        del = load_vcd_diff('2005-06-01','2005-06-30',diff_type,'pix');
        figure;
        histsemilog(del(:),50);
        if strcmpi(diff_type,'abs')
            xlabel('VCD(BRDF) - VCD(Black sky)');
        else
            xlabel('%\Delta VCD(BRDF) - VCD(Black sky)');
        end
    end

    function all_discover_validation()
        % Load the validation results using the DISCOVER data and plot the
        % correlation between BEHR and aircraft data for all four campaigns
        % and the black sky and BRDF albedo results
        
        validation_file = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/BRDF Albedo/Std-BRDF-Discover-Comparison.mat';
        V=load(validation_file);
        
        discover_validation_plots(V.D_MD,'Maryland');
        discover_validation_plots(V.D_CA,'California');
        discover_validation_plots(V.D_TX,'Texas');
        discover_validation_plots(V.D_CO,'Colorado');
        
        % We'll also make one for all four campaigns
        All.StdRet.air_no2 = cat(1, V.D_MD.StdRet.air_no2(:), V.D_CA.StdRet.air_no2(:), V.D_TX.StdRet.air_no2(:), V.D_CO.StdRet.air_no2(:));
        All.StdRet.behr_no2 = cat(1, V.D_MD.StdRet.behr_no2(:), V.D_CA.StdRet.behr_no2(:), V.D_TX.StdRet.behr_no2(:), V.D_CO.StdRet.behr_no2(:));
        All.BRDFRet.air_no2 = cat(1, V.D_MD.BRDFRet.air_no2(:), V.D_CA.BRDFRet.air_no2(:), V.D_TX.BRDFRet.air_no2(:), V.D_CO.BRDFRet.air_no2(:));
        All.BRDFRet.behr_no2 = cat(1, V.D_MD.BRDFRet.behr_no2(:), V.D_CA.BRDFRet.behr_no2(:), V.D_TX.BRDFRet.behr_no2(:), V.D_CO.BRDFRet.behr_no2(:));
        
        discover_validation_plots(All, 'All campaigns');
    end

    function discover_validation_plots(D, campaign_name)
        % Makes the individual discover-aq plots. Pass one of the
        % structures loaded in all_discover_validation.
        xcoords = [0, max(D.StdRet.air_no2(:))];
        [stdret_fit.x, stdret_fit.y, stdret_fit.legstr, stdret_fit.lineData] = calc_fit_line(D.StdRet.air_no2, D.StdRet.behr_no2, 'regression', 'rma','xcoord',xcoords);
        [brdfret_fit.x, brdfret_fit.y, brdfret_fit.legstr, brdfret_fit.lineData] = calc_fit_line(D.BRDFRet.air_no2, D.BRDFRet.behr_no2, 'regression', 'rma','xcoord',xcoords);
        
        l=gobjects(5,1);
        legstr = cell(5,1);
        figure; 
        
        l(1) = line(D.StdRet.air_no2(:), D.StdRet.behr_no2(:), 'marker', 'o', 'linewidth', 1, 'color', 'r','linestyle', 'none');
        l(2) = line(stdret_fit.x, stdret_fit.y, 'linestyle', '--', 'linewidth', 2, 'color', 'r');
        l(3) = line(D.BRDFRet.air_no2(:), D.BRDFRet.behr_no2(:), 'marker', 'd', 'linewidth', 1, 'color', 'b', 'linestyle', 'none');
        l(4) = line(brdfret_fit.x, brdfret_fit.y, 'linestyle', '--', 'linewidth', 2, 'color', 'b');
        l(5) = line(xcoords, xcoords, 'linestyle', ':', 'linewidth', 2, 'color', 'k');
        legstr{1} = 'Black-sky albedo';
        legstr{2} = stdret_fit.legstr;
        legstr{3} = 'BRDF albedo';
        legstr{4} = brdfret_fit.legstr;
        legstr{5} = '1:1';
        
        legend(l,legstr);
        set(gca,'fontsize',16);
        title(campaign_name);
        xylims;
    end

    function del = load_vcd_diff(start_date, end_date, diff_type, pix_or_grid)
        % Load each of the files for the standard and BRDF product in turn and
        % reject pixels as needed. Return either the absolute or percent
        % difference in VCDs.

        if strcmpi(pix_or_grid,'pix')
            varname = 'Data';
            catdim = 1;
        elseif strcmpi(pix_or_grid,'grid')
            varname = 'OMI';
            catdim = 3;
        else
            E.badinput('PIX_OR_GRID must be ''pix'' or ''grid''')
        end

        if ~ismember(lower(diff_type),{'abs','per'});
            E.badinput('DIFF_TYPE must be either ''abs'' or ''per''');
        end

        brdf_vcd = [];
        blacksky_vcd = [];
        dvec = datenum(start_date):datenum(end_date);
        for d=1:numel(dvec)
            fname = sprintf('OMI_BEHR_%s_%s.mat', BEHR_version, datestr(dvec(d),'yyyymmdd'));
            fprintf('Loading %s\n',fname);
            BRDF = load(fullfile(brdf_dir, fname), varname);
            BS = load(fullfile(behr_paths.behr_mat_dir,fname), varname);
            BRDF = BRDF.(varname);
            BS = BS.(varname);

            for a=1:numel(BRDF)
                % Areaweight field may not be included in Data (native pixel) structure,
                % but resetting it in the OMI structure won't hurt since we're only
                % going to check that it is 0 to see if the pixel is invalid
                BRDF(a).Areaweight = ones(size(BRDF(a).Longitude));
                BRDF(a) = omi_pixel_reject(BRDF(a),'omi',0.2,'XTrackFlags');
                BRDF(a).BEHRColumnAmountNO2Trop(BRDF(a).Areaweight==0)=nan;
                brdf_vcd = cat(catdim, brdf_vcd, BRDF(a).BEHRColumnAmountNO2Trop);
                
                BS(a).Areaweight = ones(size(BS(a).Longitude));
                BS(a) = omi_pixel_reject(BS(a),'omi',0.2,'XTrackFlags');
                BS(a).BEHRColumnAmountNO2Trop(BS(a).Areaweight==0)=nan;
                blacksky_vcd = cat(catdim, blacksky_vcd, BS(a).BEHRColumnAmountNO2Trop);
            end
        end

        if strcmpi(diff_type,'abs')
            del = brdf_vcd - blacksky_vcd;
        elseif strcmpi(diff_type,'per')
            del = reldiff(brdf_vcd, blacksky_vcd)*100;
        else
            E.notimplemented(diff_type);
        end
    end

    
end

