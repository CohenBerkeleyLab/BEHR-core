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
fprintf(fid,'\n\nMYD06:\n');
readData(Data_Chk.MYD06,fid);
fprintf(fid,'\n\nMCD43C3:\n');
readData(Data_Chk.MCD43C3,fid);

fclose(fid);

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
        years_numstr = fields{y}(i:i+3);
        fprintf(fid,'\t%s:\n',years_numstr);
        data_chk = Data.(fields{y});
        if ~iscell(data_chk)
            % Each year will contain a cell array unless all the data is
            % present, in which case, the field will be the string 'All
            % data found'. Since that's not a cell and we only care about
            % missing data, we'll skip that year.
            continue
        end
        data_by_state = struct('State',states,'Dates',{cell(size(data_chk))});
        
        for a=1:numel(states)
            A=1;
            % i will indicate were the state ends and the date range string
            % begins so that we can trim off the state.  Each entry in
            % data_chk should have the form: <state>: YYYY-MM-DD to
            % YYYY-MM-DD so the +2 accounts for the colon and space.
            i = length(states{a}) + 2;
            
            % Now go through the entries and determine which state they
            % belong to. We'll organize them into the data_by_state
            % structure so that we can print them grouped by state.
            for b=1:numel(data_chk)
                date_range = data_chk{b}; 
                if ~isempty(regexpi(date_range,states{a}))
                    data_by_state(a).Dates{A} = strtrim(date_range(i:end));
                    A = A+1;
                end
            end
            
            data_by_state(a).Dates = data_by_state(a).Dates(1:A-1);
        end
        
        % Now handle the actual printing to the file.
        for c=1:numel(data_by_state)
            fprintf(fid,'\t\t%s:\n',data_by_state(c).State);
            for d=1:numel(data_by_state(c).Dates)
                fprintf(fid,'\t\t  %s\n',data_by_state(c).Dates{d});
            end
        end
    end
end