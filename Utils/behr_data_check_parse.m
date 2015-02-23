function [  ] = behr_data_check_parse( Data_Chk )
%behr_data_check_parse Writes a text file with missing BEHR data.
%   Takes the data structure from behr_data_check and writes out a text
%   file in the current directory containing the information.  


E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that a structure was passed.
if ~isstruct(Data_Chk)
    error(E.badinput('This function only accepts a structure.'));
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN PROGRAM %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('BEHR_Data_Check.txt','w');

fprintf(fid,'OMNO2:\n');
readData(Data_Chk.OMNO2,fid);
fprintf(fid,'MYD06:\n');
readData(Data_Chk.MYD06,fid);
fprintf(fid,'MCD43C3:\n');
readData(Data_Chk.MCD43C3,fid);


end

function readData(Data,fid)
    % These are the possible states that the data can be in. Missing means
    % altogether not there, incomplete means some of the swaths aren't
    % there but at least one from each day is.
    states = {'Missing','Incomplete'};

    fields = fieldnames(Data);
    % The field names should contain the year.  Look for 4 numbers in a
    % row, and read just those. That will be used as the second level
    % header in the file.
    years = cell(size(fields));
    for y=1:numel(years)
        i = regexp(fields{y},'\d\d\d\d');
        years{i} = fields{y}(i:i+3);
    end
end