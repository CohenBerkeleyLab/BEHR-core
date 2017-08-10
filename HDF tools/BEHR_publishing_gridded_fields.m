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
        
        reprocessed_gridded_vars = {'InSituAMF', 'BEHR_R_ColumnAmountNO2Trop'};
        
        publish_only_gridded_vars = {'Latitude', 'Longitude','XTrackQualityFlags',...
            'VcdQualityFlags','BEHRQualityFlags'};
    end
    
    methods(Static = true)
        function all_vars = all_gridded_vars()
            all_vars = [BEHR_publishing_gridded_fields.cvm_gridded_vars,...
                BEHR_publishing_gridded_fields.psm_gridded_vars,...
                BEHR_publishing_gridded_fields.reprocessed_gridded_vars,... 
                BEHR_publishing_gridded_fields.publish_only_gridded_vars];
        end
    end
    
end

