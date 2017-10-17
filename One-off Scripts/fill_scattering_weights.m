% fix to make SW all nans if they're all nan or 1E-30
dates = datenum('2013-08-01'):datenum('2013-09-30');
p = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    year_d = curr_date(1:4);
    month_d = curr_date(6:7);
    day_d = curr_date(9:10);
    fprintf('%s\n',curr_date);
    fname = sprintf('OMI_BEHR_omiCloudAMF_%s%s%s.mat',year_d,month_d,day_d);
    
    load(fullfile(p,fname));
    do_i_save = false;
    for s=1:numel(Data)
        sz = size(Data(s).BEHRScatteringWeights);
        for i=1:sz(1)
            if all(Data(s).BEHRScatteringWeights(i,:) < 1e-29 | isnan(Data(s).BEHRScatteringWeights(i,:))) && ~all(isnan(Data(s).BEHRScatteringWeights(i,:)))
                Data(s).BEHRScatteringWeights(i,:) = nan;
                do_i_save = true;
            end
        end
    end
    
    if do_i_save
        fprintf('\tSaving %s\n',fullfile(p,fname));
        save(fullfile(p,fname),'Data','OMI');
    else
        fprintf('\tNo need to save %s\n',fullfile(p,fname));
    end
end