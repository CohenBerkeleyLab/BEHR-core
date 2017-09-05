classdef misc_gridding_plots
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        psm_test_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests';
        psm_newprof_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_newprofiles';
        psm_fixedfills_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_newprofiles_fixed_fills';
        psm_itest_save_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_interface_tests'
        psm_alb_unscale_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_no_albedo_scaling';
        psm_alb_scale_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_albedo_scaling';
        psm_alb_coart_dir = '/Volumes/share-sat/SAT/BEHR/PSM_Tests_coart_sea';
    end
    
    methods(Static = true)
        function make_psm_period()
            data_dir = misc_gridding_plots.psm_alb_coart_dir;
            G = GitChecker;
            G.Strict = true;
            G.addReqCommits(behr_repo_dir, '80fc4b3b'); % commit where I'd merged in the albedo changes and added the BEHR quality flags field
            G.checkState();
            
            do_overwrite = ask_yn('If files exist, overwrite them?');
            
            start_date = '2013-06-01';
            end_date = '2013-08-31';
            
            read_omno2_v_aug2012('start', start_date, 'end', end_date, 'sp_mat_dir', data_dir, 'overwrite', do_overwrite);
            BEHR_main('start', start_date, 'end', end_date, 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir, 'overwrite', do_overwrite);
            BEHR_publishing_v2('start', start_date, 'end', end_date, 'output_type', 'hdf', 'pixel_type', 'native',...
                'mat_dir', data_dir, 'save_dir', fullfile(data_dir,'native'), 'organize', false, 'overwrite', do_overwrite);
            BEHR_publishing_v2('start', start_date, 'end', start_date, 'output_type', 'hdf', 'pixel_type', 'gridded',...
                'mat_dir', data_dir, 'save_dir', fullfile(data_dir,'gridded'), 'organize', false, 'overwrite', do_overwrite);
        end
        
        function test_gridding_interface()
            load_dir = misc_gridding_plots.psm_fixedfills_dir;
            save_dir = misc_gridding_plots.psm_itest_save_dir;
            behr_grid = GlobeGrid(0.05, 'domain', 'us');
            F = dir(fullfile(load_dir, 'OMI_BEHR*.mat'));
            for a=1:numel(F)
                fprintf('Gridding BEHR data in %s\n', F(a).name);
                D = load(fullfile(load_dir, F(a).name), 'Data');
                OMI_PSM = psm_wrapper(D.Data, behr_grid, 2); %#ok<NASGU>
                save_name = strrep(F(a).name, 'BEHR', 'BEHR_PSM');
                save(fullfile(save_dir, save_name), 'OMI_PSM');
            end
        end
    end
    
end

