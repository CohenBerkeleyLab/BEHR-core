function [ behr_flags ] = behr_quality_flags( behr_amfs, behr_vis_amfs, vcd_flags, xtrack_flags, modis_ocean_flag )
%BEHR_QUALITY_FLAGS Create the BEHRQualityFlags field
%   In order to simplify things for the end user, I decided to combine the
%   NASA quality flags with some of our own into a single flag field. Like
%   VcdQualityFlags and XTrackFlags, the idea is that we will use a bit
%   array, i.e. a number where each bit that makes up its binary
%   representation has a specific meaning. 
%
%   The first (least significant) bit is reserved as an ERROR summary flag,
%   i.e. it will be == 1 if a problem with the pixel means that it should
%   not be used under normal circumstances. The idea is that end users will
%   be able to remove any pixel where this bit is set (which means the flag
%   value will be odd) and that will be sufficient filtering to remove any
%   bad pixels.
%
%   The second-sixteenth bits are to provide more information about the
%   root cause of the error, each one should be set if a certain condition
%   is met. If any of these are set, then the first bit (summary bit)
%   should also be set.
%
%   The seventeenth-thirty second bits indicate warnings, these could be
%   something that indicates that there MIGHT be a problem with the pixel,
%   but that it is usually still useable, or just noting that a certain
%   behavior occured in the retrieval for that pixel; for instance, I
%   intend to use this to mark instances where the BRDF albedo used a water
%   model, rather than the MODIS parameters.
%
%   Usage:
%
%   BEHR_FLAGS = BEHR_QUALITY_FLAGS( BEHR_AMFS, BEHR_VIS_AMFS, VCD_FLAGS,
%   XTRACK_FLAGS, MODIS_OCEAN_FLAGS )
%
%       BEHR_FLAGS: the array of flags as unsigned 32 bit integers; an
%       array the same size as BEHR_AMFS.
%
%       BEHR_AMFS: the array of BEHR total AMF values.
%
%       BEHR_VIS_AMFS: the array of BEHR visible-only AMF values.
%
%       VCD_FLAGS: the VcdQualityFlags read from NASA SP.
%
%       XTRACK_FLAGS: the XTrackQualityFlags read from NASA SP.
%
%       MODIS_OCEAN_FLAG: a logical array that is true where an ocean model
%       was used instead of the MODIS kernels and coefficients.

behr_flags = uint32(zeros(size(behr_amfs)));

%%%%%%%%%%%%%%%
% ERROR FLAGS %
%%%%%%%%%%%%%%%

% Set an error flag if the AMF has been set to the minimum value
behr_flags = set_flags(behr_flags, behr_amfs <= behr_min_amf_val() || behr_vis_amfs <= behr_min_amf_val(), 2, false);

% Set an error flag if the VcdQualityFlags field is not an even value (it's
% own quality summary flag was set)
behr_flags = set_flags(behr_flags, mod(vcd_flags, 2) ~= 0, 3, false);

% Set an error flag if the XTrackQualityFlags fields is ~= 0, i.e. it has
% been affected by the row anomaly
behr_flags = set_flags(behr_flags, xtrack_flags ~= 0, 4, false);


%%%%%%%%%%%%%%%%%
% WARNING FLAGS %
%%%%%%%%%%%%%%%%%

% Set a warning flag if we have to use an ocean model for the BRDF
behr_flags = set_flags(behr_flags, modis_ocean_flag, 9, true);

end

function flags = set_flags(flags, bool_mask, bit, warning_only)
% This subfunction should always be used to set the flags.
%
%   FLAGS: the array of BEHR quality flags to modify.
%
%   BOOL_MASK: a logical array the same size as FLAGS that is true where
%   the bit should be set in the FLAGS array.
%
%   BIT: the index of the bit to set to 1 (1-based, starting from least
%   significant)
%
%   WARNING_ONLY: a scalar logical that, if true, indicates that this bit
%   is only a warning flag. When true, it will not set the summary bit to 1
%   everywhere that BIT is also set to true. It is intended that bits 2-16
%   are error bits (i.e. WARNING_ONLY = false) and 17-32 are warning bits
%   (WARNING_ONLY = true). A warning will be issued if the bit number and
%   WARNING_ONLY do not agree, this is a check to make sure any future bits
%   follow this convention.

E = JLLErrors;
% Type and size checking 
if ~isinteger(flags)
    E.badinput('FLAGS must be an integer type')
elseif ~islogical(bool_mask) || ~isequal(size(bool_mask), size(flags))
    E.badinput('BOOL_MASK must be a logicial type the same size as FLAGS')
elseif ~isscalar(bit) || ~isnumeric(bit) || bit < 1
    E.badinput('BIT must be a positive scalar number')
elseif ~isscalar(warning_only) || ~islogical(warning_only)
    E.badinput('WARNING_ONLY must be a scalar logical');
end

if bit < 2
    E.badinput('BIT < 2: the first bit is reserved for summary flags')
elseif (bit > 16 && ~warning_only) || (bit <= 16 && warning_only)
    warning('It is expected that bits 2-16 are error bits and 17-32 are warning bits; if you are following this, be sure WARNING_ONLY reflects it');
end

flags = bitset(flags, bit, bool_mask);
if ~warning_only
    flags = bitset(flags, 1, bool_mask);
end

end