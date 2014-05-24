%Add_old_NO2
%   Adds old NO2 tropospheric columns to the new SP data sets. Solely for
%   testing purposes; this is intended to allow the revised version
%   BEHR_main.m to calculate the BEHR NO2 column with the same OMI NO2
%   column as Ashley used.  This will allow comparison to see where the
%   error is coming from - different OMI NO2 columns or my values for
%   GLOBETerpres, MODISAlbedo, etc.
%
%   Josh Laughner <joshlaugh5@gmail.com> 25 Apr 2014

%Dates in yyyy/mm/dd format
date_start = '2011/04/02';
date_end = '2011/09/30';

%The directories the files can be found in
old_sp_dir = '/Volumes/share/GROUP/SAT/BEHR/SP_files';
new_sp_dir = '/Volumes/share/GROUP/SAT/BEHR/Test_SP_files';

%The file name up to the date
old_prefix = 'OMI_SP_';
new_prefix = 'OMI_SP_';

dates = datenum(date_start):datenum(date_end);
ndates = numel(dates);

for d=1:ndates
    date = datestr(dates(d),29);
    fprintf('\n %s: \n',date);
    year = date(1:4);
    month = date(6:7);
    day = date(9:10);
    
    load(fullfile(old_sp_dir,[old_prefix,year,month,day,'.mat']));
    oldData = Data; clear Data
    
    load(fullfile(new_sp_dir,[new_prefix,year,month,day,'.mat']));
    newData = Data; clear Data
    
    compare = Compare_Data_Fields(oldData,newData,'ColumnAmountNO2Trop'); %Use this simply to find the pixels that match, the field doesn't matter
    matched_pix = {};
    for a=1:length(compare)
        matched_pix = [matched_pix; compare(a).Array(2:end,:)];
    end
    % First, copy all the fields from newData in the indices that match up
    % with old pixels
    
    clear appendedData
    new_indicies = unique(cell2mat(matched_pix(2:end,3)));
    index_diff = min(new_indicies)-1; %This will ensure that even if the new indices are e.g. [2, 3, 4] the appended structure will have indices [1, 2, 3]
    
    fprintf('   Copying new fields\n');
    for a=1:length(new_indicies)
        appendedData(a) = newData(new_indicies(a));
    end
    
    % Then add the specific fields from the old data we want to
    
    fprintf('   Copying old NO2 column and AMFTrop\n');
    for a=1:size(matched_pix,1)
        A_ind = matched_pix{a,3} - index_diff; %Index for the appended data structure
        O_ind = matched_pix{a,1}; %Index for the old data structure
        A_pix = matched_pix{a,4}; %Pixel for the appended data structure (same as the new data. str.)
        O_pix = matched_pix{a,2}; %Pixel for the old data structure
        
        appendedData(A_ind).OldNO2ColumnTrop(A_pix) = oldData(O_ind).ColumnAmountNO2Trop(O_pix);
        appendedData(A_ind).OldAMFTrop(A_pix) = oldData(O_ind).AMFTrop(O_pix);
        appendedData(A_ind).OldGLOBETerpres(A_pix) = oldData(O_ind).GLOBETerpres(O_pix);
        appendedData(A_ind).OldMODISCloud(A_pix) = oldData(O_ind).MODISCloud(O_pix);
    end
    
    save(fullfile(new_sp_dir,[new_prefix,'appendedGLOBE_',year,month,day,'.mat']),'appendedData')
    
end
