function [ wrf_path ] = find_wrf_path( profile_mode, this_date )
%FIND_WRF_PATH Returns the path to WRF profiles for the given date
%   WRF_PATH = FIND_WRF_PATH( PROFILE_MODE, THIS_DATE ) Returns the proper
%   path to the WRF-Chem profiles as WRF_PATH. Given PROFILE_MODE =
%   'monthly', will just return behr_paths.wrf_monthly_profiles. Given
%   PROFILE_MODE = 'daily', it will search all paths defined in the cell
%   array behr_paths.wrf_profiles for the year and month subdirectories
%   that match THIS_DATE. If it finds one, it will return it. Otherwise,
%   an error is thrown. It does not verify that the required wrfout files
%   are actually present. THIS_DATE must be either a date number or date
%   string implicitly understood by Matlab.

E = JLLErrors;

validate_date(this_date);

if strcmpi(profile_mode, 'monthly')
    wrf_path = behr_paths.wrf_monthly_profiles;
    return
elseif strcmpi(profile_mode, 'daily')
    yr_str = datestr(this_date, 'yyyy');
    mn_str = datestr(this_date, 'mm');
    wrf_dirs = behr_paths.wrf_profiles;
    for a=1:numel(wrf_dirs)
        if ~exist(wrf_dirs{a}, 'dir')
            E.dir_dne('The root WRF-Chem profile directory %s does not exist;\n  Is the file server mounted?\n  Is the directory defined correctly in behr_paths.m?', wrf_dirs{a});
        end
        
        wrf_dir = fullfile(wrf_dirs, yr_str, mn_str);
        if exist(wrf_dir, 'dir')
            return
        end
    end
    
    % If we've gotten here, we haven't found the directory
    E.dir_dne('No WRF-Chem daily output directory exists for %s', datestr(this_date, 'mmm yyyy'));
else
    E.badinput('No paths defined for PROFILE_MODE = ''%s''', profile_mode);
end

end

