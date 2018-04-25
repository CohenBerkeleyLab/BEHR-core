function [ success ] = main_priori_tests( data, DEBUG_LEVEL, fid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('DEBUG_LEVEL', 'var')
    DEBUG_LEVEL = 2;
end

if ~exist('fid', 'var')
    % An fid of 1 will make fprint print to the command window as if no fid
    % was given
    fid = 1;
end

myname = mfilename();
success = check_quality_flags(data);

    function passfail = check_quality_flags(data)
        % Verify that:
        %   * The error summary bit is true everywhere that the pixel
        %     absolutely should not be used.
        %   * The quality summary bit is true everywhere the error bit is
        %   * The quality summary bit is also true everywhere that we
        %     declare the to-ground column quality will be low
        % See each part of the test for exactly what meets these criteria.
        
        if DEBUG_LEVEL > 0
            fprintf(fid, 'Checking validity of quality flags...\n');
        end
        passfail = true;
        
        for a=1:numel(data)
            quality_bit = bitand(data(a).BEHRQualityFlags, 1) > 0;
            error_bit = bitand(data(a).BEHRQualityFlags, 2) > 0;
            
            % First verify that the error bit is correct. As of v3.0B, this
            % should be set anywhere that either of the BEHR AMFs is the
            % minimum value or a NaN, anywhere that the VcdQualityFlags
            % indicate and error in NASA processing, or anywhere that
            % XTrackQualityFlags indicates that we're in the row anomaly. If
            % this changes in the future, be sure to update this test to match
            % what is expected.
            error_check = data(a).BEHRAMFTrop <= behr_min_amf_val | isnan(data(a).BEHRAMFTrop) | data(a).BEHRAMFTropVisOnly <= behr_min_amf_val | isnan(data(a).BEHRAMFTropVisOnly)...
                | mod(data(a).VcdQualityFlags, 2) ~= 0 | data(a).XTrackQualityFlags > 0;
            
            error_mismatch = xor(error_bit(:), error_check(:));
            
            
            % Next verify that the quality bit is true everywhere the
            % error bit is. This test should be the same for any future
            % versions of BEHR.
            qual_error_mismatch = error_bit(:) & ~quality_bit(:);
            
            % Finally, verify that in addition to the error bit, the quality
            % fit flags low quality pixels. As of v3.0B, this should include
            % pixels with CloudFraction > 0.2 and low quality MODIS data
            quality_check = error_bit | data(a).CloudFraction > 0.2 | data(a).MODISAlbedoQuality >= 2.5 | data(a).MODISAlbedoFillFlag;
            quality_mismatch = xor(quality_bit(:), quality_check(:));
            
            % Now summarize for this orbit and print out any necessary
            % messages (true == pass)
            this_passfail = ~any(error_mismatch) && ~any(qual_error_mismatch) && ~any(quality_mismatch);
            if ~this_passfail
                % If any orbit fails, fail the overall test, but we can't
                % just do passfail = this_passfail, since then the overall
                % test will just represent the last orbit.
                passfail = false;
            end
            
            if DEBUG_LEVEL > 0
                fprintf(fid, '  Swath %d: %s\n', a, passfail_str(this_passfail));
            end
            if DEBUG_LEVEL > 1
                % Specific reasons
                if any(error_mismatch)
                    fprintf(fid, '    %1$d pixels'' error bits do not agree with what we expect. If you have changed the definition of the error bit in behr_quality_flags but not %2$s, you may need to update %2$s\n', sum(error_mismatch), myname);
                end
                
                if any(qual_error_mismatch)
                    fprintf(fid, '    %d pixels have the error bit but not the quality bit set.\n', sum(qual_error_mismatch));
                end
                
                if any(quality_mismatch)
                    fprintf(fid, '    %1$d pixels'' quality bits do not agree with what we expect. If you have changed the definition of the quality bit in behr_quality_flags but not %2$s, you may need to update %2$s\n', sum(quality_mismatch), myname);
                end
            end
            
        end
    end

end

function s = passfail_str(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end