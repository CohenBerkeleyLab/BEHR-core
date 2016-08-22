function [  ] = check_fit_goodness( )
%CHECK_FIT_GOODNESS Prints out a list of goodness of fit metrics for fits
%   Opens dialogue boxes to choose the fit files and the fit metrics to
%   print. 

E = JLLErrors;
[data_files, data_path] = uigetfile('*SimpleFits*.mat','Select the fits files to print goodness-of-fit metrics for','multiselect','on');
if isnumeric(data_files) && data_files == 0
    E.userCancel;
elseif ischar(data_files)
    data_files = {data_files};
end

F = load(fullfile(data_path, data_files{1}));
stats = fieldnames(F.f_hyfast.stats);

stat_ind = listdlg('ListString',stats);
stat = stats{stat_ind};

for a=1:numel(data_files)
    F = load(fullfile(data_path, data_files{a}));
    fns = fieldnames(F);
    city_name = regexp(data_files{a},'(Atlanta|Birmingham|Montgomery)', 'match', 'once');
    wind_spd = regexp(data_files{a}, '\dpt\d', 'match', 'once');
    for b=1:numel(fns)
        s = F.(fns{b}).stats.(stat);
        if isscalar(s)
            valstr = sprintf('%s = %.4f', stat, s);
        elseif numel(s) == 5
            valstr = sprintf('%1$s a = %2$g, %1$s x0 = %3$g, %1$s mux = %4$g, %1$s sigma_x = %5$g, %1$s B = %6$g', stat, s(1), s(2), s(3), s(4), s(5));
        else
            E.notimplemented('stat not a 5 element vector or scalar')
        end
        fprintf('%s %s %s: %s\n', city_name, wind_spd, fns{b}, valstr);
    end
end

end

