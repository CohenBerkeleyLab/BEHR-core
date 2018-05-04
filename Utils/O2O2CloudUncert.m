classdef O2O2CloudUncert < handle
    %O2O2CloudUncert Stores cloud pressure uncertainty from Acarreta et al. 2008
    %   O2 = O2O2CloudUncert() will automatically read in the text file
    %   containing the GraphClick-extracted fig. 3 from Acarreta et al.
    %   (doi:10.1029/2003JD003915) and create a scattered interpolant
    %   (O2.interpolant) representing that data.
    
    properties
        error_vals;
        cld_frac;
        cld_pres;
        interpolant;
    end
    
    methods
        function obj = O2O2CloudUncert()
            %O2O2CloudUncert Construct an instance of this class
            obj.parse_graphclick_file;
        end
        
        function parse_graphclick_file(obj)
            %PARSE_GRAPHCLICK_FILE Parses the .txt file created by GraphClick
            %   I extracted the data from Acarreta et al.'s Fig. 3 using
            %   GraphClick (http://www.arizona-software.ch/graphclick/) and
            %   stored it as the text file (whose path is returned by
            %   O2O2CloudUncert.figure_text_file). This function reads that
            %   file and generates the necessary interpolant.
            fid = fopen(obj.figure_text_file);
            tline = fgetl(fid);
            i_val = 1;
            curr_error = nan;
            
            n_val = 3000; % there's about 3000 values in the text file as of 23 Apr 2018
            these_error_vals = nan(n_val,1);
            these_cld_frac = nan(n_val,1);
            these_cld_pres = nan(n_val,1);
            while ischar(tline)
                % File will have the format:
                %   <isopleth>:
                %   x       y
                %   <data>
                % where the error in cloud pressure is the given isopleth
                % and the X and Y data give the cloud fraction and cloud
                % pressure (respectively) corresponding to that error. We
                % need to convert these into a scattered interpolant so
                % that we can query arbitrary cloud pressure/cloud
                % fractions.
                if regcmp(tline, ':')
                    curr_error = str2double(regexp(tline, '\d+', 'match', 'once'));
                elseif ~regcmp(tline, 'x\s+y') && ~isempty(tline)
                    if isnan(curr_error)
                        error('O2O2CloudUncert:parse_error', 'Current error is a NaN while trying to read coordinate data.');
                    end
                    vals = cellfun(@str2double, strsplit(tline));
                    these_error_vals(i_val) = curr_error;
                    these_cld_frac(i_val) = vals(1);
                    these_cld_pres(i_val) = vals(2);
                    i_val = i_val + 1;
                end
                
                tline = fgetl(fid);
            end
            fclose(fid);
            
            obj.error_vals = these_error_vals(1:(i_val-1));
            obj.cld_frac = these_cld_frac(1:(i_val-1));
            obj.cld_pres = these_cld_pres(1:(i_val-1));
            obj.interpolant = scatteredInterpolant(obj.cld_frac, obj.cld_pres, obj.error_vals);
        end
    end
    
    methods(Static)
        function value = figure_text_file()
            % FIGURE_TEXT_FILE Return the path to the GraphClick-generated text file that this reads.
            mydir = fileparts(mfilename('fullpath'));
            value = fullfile(mydir, 'Constants', 'ErrorAnalysis', 'AcarrataCldPresUncert.txt');
        end
        
    end
end

