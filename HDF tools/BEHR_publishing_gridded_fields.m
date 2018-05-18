classdef BEHR_publishing_gridded_fields
    %BEHR_PUBLISHING_GRIDDED_FIELDS Representation of which data fields should be gridded in the BEHR algorithm
    %   Starting from version 3 of BEHR, we have the capability to grid our
    %   data with either the constant value or parabolic spline method (CVM
    %   and PSM respectively) from Gerrit Kuhlmann's omi Python package. In
    %   keeping with the programming paradigm that all data should have a
    %   single, authoritative representation within a program, this class
    %   indicates which data fields should be gridded using each algorithm.
    %
    %   Properties (all cell arrays of strings):
    %       cvm_gridded_vars: standard data fields that should be gridded
    %       using the CVM method.
    %
    %       psm_gridded_vars: standard data fields that can be gridded
    %       using the PSM method (they may also be gridded with CVM, if
    %       desirable).
    %
    %       reprocessed_psm_vars; reprocessed_cvm_vars: analogously to the
    %       first two, theses are variables created during the reprocessing
    %       of BEHR data with in situ measured profiles.
    %
    %       flag_vars: this are variables that should be added using a
    %       bitwise OR operation rather than the usual weighted average.
    %       This assumes flag variables are usually a bit array, where each
    %       bit that makes up a value in binary has some meaning.
    %
    %       publish_only_gridded_vars: variables that do not get passed to
    %       the omi gridding package, but are still present in the gridded
    %       data and should be included in any published files. Usually
    %       these define the grid (i.e. latitude/longitude).
    %
    %   Methods:
    %       all_gridded_vars(): returns a cell array with all the gridded
    %       variables included.
    %
    %       all_psm_vars(): returns a cell array with all the fields that
    %       can be gridded by PSM (so the standard + reprocessed PSM
    %       fields).
    %
    %       psm_weight_vars(): returns a cell array with the weight field
    %       names corresponding to all of the PSM fields. The gridding
    %       interface allows for each PSM gridded field to have different
    %       weights, so individual weight fields are necessary.
    %
    %       make_wts_field(): given a field name, returns the corresponding
    %       weights field name, assuming it is a PSM variable.
    
    properties(Constant = true)
        cvm_gridded_vars = {'AmfTrop', 'BEHRAMFTrop','BEHRAMFTropVisOnly',...
            'CloudFraction', 'CloudRadianceFraction',  'BEHRSurfacePressure',...
            'MODISAlbedo', 'MODISCloud', 'Row','BEHRTropopausePressure'};
       
        psm_gridded_vars = {'BEHRColumnAmountNO2Trop',...
            'BEHRColumnAmountNO2TropVisOnly',...
            'ColumnAmountNO2Trop'};
        
        % Currently I have no mechanism to really handle these. If, in the
        % future, I reprocess data around the DISCOVER campaigns with
        % version 3, then these will be used to set which variables created
        % using the in situ measured profiles are gridded with PSM and
        % which with CVM.
        reprocessed_psm_vars = {'BEHR_R_ColumnAmountNO2Trop'};
        reprocessed_cvm_vars = {'InSituAMF'};
        
        % These are variables that will be gridded by CVM using a bitwise
        % OR operation, rather than a weighted average.
        flag_vars = {'XTrackQualityFlags', 'VcdQualityFlags', 'BEHRQualityFlags'};
        
        % These are variables that will be included in the published data,
        % but aren't gridded normally. Usually, as in the case of latitude
        % and longitude, they define the grid itself.
        publish_only_gridded_vars = {'Latitude', 'Longitude'};

        % These are fields that go into each swath's attributes, so they need to be
        % copied by psm_wrapper into the OMI structure
        swath_attr_vars = {'Date', 'OMNO2File', 'OMPIXCORFile', 'MODISCloudFiles', 'MODISAlbedoFile',...
                           'BEHRWRFFile', 'BEHRProfileMode', 'BEHRRegion', 'BEHRWRFPressureMode',...
                           'BEHRWRFTemperatureMode', 'GitHead_Core_Read', 'GitHead_BEHRUtils_Read',...
                           'GitHead_GenUtils_Read', 'GitHead_Core_Main', 'GitHead_BEHRUtils_Main',...
                           'GitHead_GenUtils_Main', 'GitHead_PSM_Main', 'GitHead_MatPyInt_Main',...
                           'GitHead_WRFUtils_Main'};
    end
    
    methods(Static = true)
        function all_vars = all_gridded_vars(subset)
            if ~exist('subset', 'var')
                subset = 'all';
            end
            
            all_vars = [BEHR_publishing_gridded_fields.cvm_gridded_vars,...
                BEHR_publishing_gridded_fields.psm_gridded_vars,...
                BEHR_publishing_gridded_fields.psm_weight_vars(subset),...
                BEHR_publishing_gridded_fields.cvm_weight_vars,...
                BEHR_publishing_gridded_fields.reprocessed_psm_vars,... 
                BEHR_publishing_gridded_fields.reprocessed_cvm_vars,...
                BEHR_publishing_gridded_fields.flag_vars,...
                BEHR_publishing_gridded_fields.publish_only_gridded_vars];
        end
        
        function all_psm_vars = all_psm_vars()
            all_psm_vars = [BEHR_publishing_gridded_fields.psm_gridded_vars,...
                BEHR_publishing_gridded_fields.reprocessed_psm_vars];
        end
        
        function all_cvm_vars = all_cvm_vars()
            all_cvm_vars = [BEHR_publishing_gridded_fields.cvm_gridded_vars,...
                BEHR_publishing_gridded_fields.reprocessed_cvm_vars];
        end
        
        function weight_vars = psm_weight_vars(subset)
            if ~exist('subset', 'var')
                subset = 'all';
            end

            if strcmpi(subset, 'std')
                weight_vars = BEHR_publishing_gridded_fields.psm_gridded_vars;
            elseif strcmpi(subset, 'insitu')
                weight_vars = BEHR_publishing_gridded_fields.reprocessed_psm_vars;
            else
                weight_vars = BEHR_publishing_gridded_fields.all_psm_vars;
            end

            for a=1:numel(weight_vars)
                weight_vars{a} = BEHR_publishing_gridded_fields.make_wts_field(weight_vars{a});
            end
        end

        function weight_vars = cvm_weight_vars()
            weight_vars = {'Areaweight'};
        end
        
        function fn = make_wts_field(fn)
            fn = sprintf('%sWeights', fn);
        end
    end
    
end

