function [ att_val ] = hdfreadatt( filename, dsetname, attname )
%HDFREADATT Read an attribute value from a given dataset
%   ATT_VAL = HDFREADATT(FILENAME, DSETNAME, ATTNAME) Read the attribute
%   with the name ATTNAME from the SDS DSETNAME in file FILENAME. DSETNAME
%   should be a / separated path of Vgroups ending in the dataset name.
%   E.g. if the dataset 'Longitude' is in the 'Geolocation Fields' Vgroup,
%   which in turn is under the 'mod06' Vgroup, DSETNAME should be
%   '/mod06/Geolocation Fields/Longitude'. The leading / is not necessary.
%
%   Why "h5readatt" is a built-in function and this isn't is beyond me.

E = JLLErrors;
if ~ischar(filename)
    E.badinput('FILENAME must be a string')
elseif ~exist(filename, 'file')
    E.badinput('%s does not exist', filename)
end

if ~ischar(dsetname)
    E.badinput('DSETNAME must be a string')
end

if ~ischar(attname)
    E.badinput('ATTNAME must be a string')
end

%%%% MAIN FUNCTION %%%%

dset_path = strsplit(dsetname, '/');
dset_path(iscellcontents(dset_path, 'isempty')) = []; % handles the empty string caused by a leading /

hinfo = hdfinfo(filename);

% Assume all but the last part of the path are Vgroups, and the last one is
% a dataset
for a=1:numel(dset_path) - 1
    vgroup_names = {hinfo.Vgroup.Name};
    dind = strcmp(vgroup_names, dset_path{a});
    if sum(dind) ~= 1
        path_so_far = strjoin(dset_path(1:a), '/');
        E.callError('hdf_vgroup', 'Cannot find %s in %s', path_so_far, filename);
    end
    hinfo = hinfo.Vgroup(dind);
end

dset_names = {hinfo.SDS.Name};
dind = strcmp(dset_names, dset_path{end});
if sum(dind) ~= 1
    E.callError('hdf_sds', 'Cannot find %s in %s', dsetname, filename);
end

dset = hinfo.SDS(dind);

dset_atts = {dset.Attributes.Name};
aind = strcmp(dset_atts, attname);
if sum(aind) ~= 1
    E.callError('hdf_attribute', 'Cannot find attribute %s in SDS %s', attname, dsetname);
end

att_val = dset.Attributes(aind).Value;

end

