classdef BEHR_publishing_gridded_fields
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        cvm_gridded_vars = {'AmfTrop', 'BEHRAMFTrop','BEHRAMFTropVisOnly',...
            'CloudFraction', 'CloudRadianceFraction',  'GLOBETerpres',...
            'MODISAlbedo', 'MODISCloud', 'Row'};
        
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
    end
    
    methods(Static = true)
        function all_vars = all_gridded_vars()
            all_vars = [BEHR_publishing_gridded_fields.cvm_gridded_vars,...
                BEHR_publishing_gridded_fields.psm_gridded_vars,...
                BEHR_publishing_gridded_fields.psm_weight_vars,...
                BEHR_publishing_gridded_fields.reprocessed_psm_vars,... 
                BEHR_publishing_gridded_fields.reprocessed_cvm_vars,...
                BEHR_publishing_gridded_fields.flag_vars,...
                BEHR_publishing_gridded_fields.publish_only_gridded_vars];
        end
        
        function all_psm_vars = all_psm_vars()
            all_psm_vars = [BEHR_publishing_gridded_fields.psm_gridded_vars,...
                BEHR_publishing_gridded_fields.reprocessed_psm_vars];
        end
        
        function weight_vars = psm_weight_vars()
            weight_vars = BEHR_publishing_gridded_fields.all_psm_vars;
            for a=1:numel(weight_vars)
                weight_vars{a} = BEHR_publishing_gridded_fields.make_wts_field(weight_vars{a});
            end
        end
        
        function fn = make_wts_field(fn)
            fn = sprintf('%sWeights', fn);
        end
    end
    
end

