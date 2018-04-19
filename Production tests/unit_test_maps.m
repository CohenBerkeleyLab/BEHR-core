function  unit_test_maps(ut_base_dir, ut_new_dir, fields, varargin)
%UNIT_TEST_MAPS Make maps of the differences between two unit tests
%   UNIT_TEST_MAPS() will interactively request all the necessary options
%
%   UNIT_TEST_MAPS( UT_BASE_DIR, UT_NEW_DIR, FIELDS ) will read the
%   OMI_BEHR .mat files from UT_BASE_DIR and UT_NEW_DIR and plot
%   differences (relative by default) for each field specified by FIELDS in
%   the Data structure.
%
%   Parameter arguments:
%
%       'diff_type' - must be the string 'rel' (default, relative percent
%       differences), 'abs' (absolute difference), 'nan' (plot the
%       difference in whether the value of the fields is a NaN), 'base'
%       (plot the base value, no difference), 'new' (plot the new value,
%       no difference), 'basenan', or 'newnan' (plot whether the base or
%       new value is a NaN, respectively).
%
%           This can also be given as a structure, where each
%       field name is one of the fields to be plotted and the field value
%       is one of the allowed values of diff_type. In this form, each field
%       can have it's own difference type, but the structure must include
%       every field.
%
%       'structure' - must be the strings 'Data' (default) or 'OMI',
%       controls which structure in the .mat file is plotted.
%
%       'close' - logical, default is true, which causes each day's plots
%       to close when you go onto the next day. False will keep each day's
%       plots open.
%
%       'clim' - specify a color limit to be passed to caxis() for each
%       plot. The default is set using calc_plot_limits( del, 'pow10',
%       'diff' ), where del is the vector of differences.
%
%       'mode_3d' - how to deal with 3D fields. Options are: 'avg', meaning
%       that differences are averaged along the vertical dimension before
%       plotting. Other options are 'absavg', which takes the absolute
%       value of the differences before averaging; 'sum', which adds up the
%       differences; and 'abssum', which adds up the absolute values of the
%       differences. For diff_type == 'rel' or 'abs', 'avg' is the default;
%       for diff_type == 'nan', 'abssum' is the default.
E = JLLErrors;

p = inputParser;
p.addParameter('diff_type', 'rel');
p.addParameter('structure', 'Data');
p.addParameter('close', true);
p.addParameter('clim', []);
p.addParameter('mode_3d','');

if nargin == 0
    [ut_base_dir, ut_new_dir, fields, varargin] = get_settings_interactive();
end

p.parse(varargin{:});
pout = p.Results;

diff_type = pout.diff_type;
structure = pout.structure;
do_close = pout.close;
clim = pout.clim;
mode_3d = pout.mode_3d;

if ~exist(ut_base_dir, 'dir')
    E.badinput('UT_BASE_DIR must be a valid directory')
end
if ~exist(ut_new_dir, 'dir')
    E.badinput('UT_NEW_DIR must be a valid directory')
end

if ischar(fields)
    fields = {fields};
elseif ~iscellstr(fields)
    E.badinput('FIELDS must be a char or cell array of chars');
end

allowed_diff_types = {'rel', 'abs', 'nan', 'base', 'new', 'basenan', 'newnan'};
if ischar(diff_type)
    if ~ismember(diff_type, allowed_diff_types)
        E.badinput('If given as a char, ''diff_type'' must be one of: %s', strjoin(allowed_diff_types, ', '))
    end
    
    diff_type = make_empty_struct_from_cell(fields, diff_type);
elseif isstruct(diff_type)
    missing_fields = fields(~isfield(diff_type, fields));
    if ~isempty(missing_fields)
        E.badinput('If given as a struct, ''diff_type'' must have every value in FIELDS as a field. The following fields are missing: %s', strjoin(missing_fields, ', '));
    elseif any(~structfun(@ischar, diff_type)) || any(~structfun(@(x) ismember(x, allowed_diff_types), diff_type))
        E.badinput('One or more of the fields in ''diff_type'' is not one of the allowed values: %s', strjoin(allowed_diff_types, ', '));
    end
end

allowed_structures = {'Data', 'OMI'};
if ~ischar(structure) || ~ismember(structure, allowed_structures)
    E.badinput('''structure'' must be one of: %s', strjoin(allowed_structures, ', '));
end

if ~islogical(do_close) || ~isscalar(do_close)
    E.badinput('''close'' must be a scalar logical')
end

if ~isempty(clim) && (~isnumeric(clim) || numel(clim) ~= 2)
    E.badinput('''clim'' must be a 2-element numeric vector, if given')
end

% We have different 3D operation modes depending on what difference type
% we're taking
if isempty(mode_3d)
    if strcmpi(diff_type, 'nan')
        mode_3d = 'abssum';
    else
        mode_3d = 'avg';
    end
end

allowed_mode_3ds = {'avg', 'absavg', 'sum', 'abssum'};
if ~ischar(mode_3d) || ~ismember(mode_3d, allowed_mode_3ds)
    E.badinput('''mode_3d'' must be one of: %s', strjoin(allowed_mode_3ds, ', '));
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

% First get the list of files in each directory for the same date
[Fbase, Fnew] = list_common_behr_files(ut_base_dir, ut_new_dir);

% We can assume the same number of files will be selected
for i_file = 1:numel(Fbase)
    figs = gobjects(size(fields));
    for i_field = 1:numel(fields)
        this_diff_type = diff_type.(fields{i_field});
        [del, lon, lat] = diff_files(Fbase(i_file).name, Fnew(i_file).name, fields{i_field}, structure, this_diff_type, mode_3d);
        figs(i_field) = figure;
        sp_size = square_subplot_dims(numel(del), 'exact');
        for i_del = 1:numel(del)
            subplot(sp_size(1), sp_size(2), i_del);
            try
                pcolor(lon{i_del}, lat{i_del}, del{i_del});
                shading flat;
                cb=colorbar;
                cb.Label.String = cb_label(fields{i_field}, this_diff_type);
                state_outlines('k');
                
                if isempty(clim)
                    fig_clim = calc_plot_limits( del{i_del}, 'pow10', diff_limit_type(this_diff_type) );
                else
                    fig_clim = clim;
                end
                caxis(fig_clim);
                title(datestr(date_from_behr_filenames(Fbase(i_file).name)))
            catch err
                if strcmpi(err.identifier, 'MATLAB:pcolor:NonMatrixColorInput') && isvector(del{i_del})
                    fprintf('Skipping orbit %d because pcolor can''t handle vector inputs\n', i_del);
                else
                    rethrow(err);
                end
            end
        end
        colormap(diff_colormap(this_diff_type));
    end
    if do_close
        tilefigs;
        input('Press ENTER to continue', 's');
        close(figs);
    end
end

end

function [ut_base_dir, ut_new_dir, fields, opts_cell] = get_settings_interactive()
ut_base_dir = ask_file('Choose the base directory with OMI_BEHR .mat files', 'dir');
ut_new_dir = ask_file('Choose the new directory with OMI_BEHR .mat files', 'dir');

opts.structure = ask_multichoice('Which structure to use?', {'Data', 'OMI'}, 'list', true, 'default', 'Data');

% load the first base file to get the list of available fields. Only allow
% the user to choose numeric fields.
F = dirff(fullfile(ut_base_dir, 'OMI_BEHR*.mat'));
tmp = load(F(1).name, opts.structure);
S = tmp.(opts.structure);
fns = fieldnames(S);
keep_fn = true(size(fns));
for i_fn = 1:numel(fns)
    keep_fn(i_fn) = isnumeric(S(1).(fns{i_fn}));
end

fields = ask_multiselect('Which field(s) to plot?', fns(keep_fn));

opts.diff_type = ask_multichoice('Which difference type to plot?', {'rel','abs','nan'}, 'list', true);
opts.close = ask_yn('Close each day''s figure before moving onto the next?');
if strcmpi(opts.diff_type, 'nan')
    default_3d = 'abssum';
else
    default_3d = 'avg';
end
opts.mode_3d = ask_multichoice('How to plot 3D fields?', {'avg', 'absavg', 'sum', 'abssum'}, 'list', true, 'default', default_3d);

opts_cell = struct2cell2(opts);

end


function [Fbase, Fnew] = list_common_behr_files(ut_base_dir, ut_new_dir)
Fbase = dirff(fullfile(ut_base_dir, 'OMI_BEHR*.mat'));
Fnew = dirff(fullfile(ut_new_dir, 'OMI_BEHR*.mat'));
base_dates = date_from_behr_filenames(Fbase);
new_dates = date_from_behr_filenames(Fnew);
[perm_base, perm_new] = find_common_elements(base_dates, new_dates);
Fbase = Fbase(perm_base);
Fnew = Fnew(perm_new);
end


function [del, lon, lat] = diff_files(base_name, new_name, field, structure, diff_type, mode_3d)
E = JLLErrors;

switch lower(diff_type)
    case 'rel'
        del_fxn = @(new, base) reldiff(new,base)*100;
    case 'abs'
        del_fxn = @(new, base) new - base;
    case 'nan'
        del_fxn = @(new, base) double(isnan(new)) - double(isnan(base));
    case 'base'
        del_fxn = @(new, base) base;
    case 'new'
        del_fxn = @(new, base) new;
    case 'basenan'
        del_fxn = @(new, base) double(isnan(base));
    case 'newnan'
        del_fxn = @(new, base) double(isnan(new));
    otherwise
        E.notimplemented('No difference function defined for diff_type == "%s"', diff_type);
end

switch lower(mode_3d)
    case 'avg'
        collapse_3d = @(x) squeeze(nanmean(x,1));
    case 'absavg'
        collapse_3d = @(x) squeeze(nanmean(abs(x),1));
    case 'sum'
        collapse_3d = @(x) squeeze(nansum2(x,1));
    case 'abssum'
        collapse_3d = @(x) squeeze(nansum2(abs(x),1));
    otherwise
        E.notimplemented('No method for reducing 3D arrays defined for mode_3d == "%s"', mode_3d);
end

Base = load(base_name, structure);
Base = Base.(structure);
New = load(new_name, structure);
New = New.(structure);

if ~isequal(size(New), size(Base))
    E.notimplemented('New and Base structures are different sizes');
end

init_cell = cell(size(Base));
del = init_cell;
lon = init_cell;
lat = init_cell;

for i = 1:numel(Base)
    del{i} = del_fxn(New(i).(field), Base(i).(field));
    
    if ~ismatrix(del{i})
        del{i} = collapse_3d(del{i});
    end
    
    if ~isequaln(New(i).Longitude, Base(i).Longitude) || ~isequaln(New(i).Latitude, Base(i).Latitude)
        E.notimplemented('New and Base have different lat/lon coordinates')
    end
    
    lon{i} = New(i).Longitude;
    lat{i} = New(i).Latitude;
end

end

function label = cb_label(field, diff_type)
E = JLLErrors;
switch lower(diff_type)
    case 'rel'
        label = sprintf('%%\\Delta %s', field);
    case 'abs'
        label = sprintf('\\Delta %s', field);
    case 'nan'
        label = sprintf('isnan(New) - isnan(Old) (%s)', field);
    case 'base'
        label = sprintf('Base %s', field);
    case 'new'
        label = sprintf('New %s', field);
    case 'basenan'
        label = sprintf('isnan(Base) (%s)', field);
    case 'newnan'
        label = sprintf('isnan(New) (%s)', field);
    otherwise
        E.notimplemented('No label defined for diff_type == "%s"', diff_type);
end
end

function lim_type = diff_limit_type(diff_type)
if regcmp(diff_type, '(base|new)')
    lim_type = 'zero';
else
    lim_type = 'diff';
end
end

function cmap = diff_colormap(diff_type)
if regcmp(diff_type, '(base|new)')
    cmap = parula;
else
    cmap = blue_red_cmap;
end
end