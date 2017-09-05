function [ output_args ] = draw_modis_cloud_swaths( year_in, day_of_year )
%draw_modis_cloud_swaths Plot the outline of MYD06 datasets

DEBUG_LEVEL = 1;
E = JLLErrors;

p = '/Volumes/share-sat/SAT/MODIS/MYD06_L2/';
if isnumeric(year_in);
    year_in = num2str(year_in);
end
if isnumeric(day_of_year)
    day_of_year = sprintf('%03d',day_of_year);
elseif ischar(day_of_year)
    day_of_year = sprintf('%03s',day_of_year);
else
    E.badinput('day_of_year must be either a string or integer')
end

F = dir(fullfile(p,year_in,sprintf('MYD06_L2.A%s%s*',year_in,day_of_year)));

% Initialize a map of the US
figure;
m_proj('Mercator','long',[-130 -60],'lat',[20 55]);
m_states('k');
m_grid('linestyle','none');

for a=1:numel(F)
    if DEBUG_LEVEL > 0
        if a > 1
            erasestr = repmat('\b',1,length(str));
            fprintf(erasestr);
        end
        str = sprintf('Now plotting file %d of %d',a,numel(F));
        fprintf(str);
        if a==numel(F)
            fprintf('\n');
        end
    end
    mlat = hdfread(fullfile(p,year_in,F(a).name),'Latitude');
    mlon = hdfread(fullfile(p,year_in,F(a).name),'Longitude');
    x = zeros(1,5); y = zeros(1,5);
    x(1) = double(mlon(1,1)); y(1) = double(mlat(1,1));
    x(2) = double(mlon(1,end)); y(2) = double(mlat(1,end));
    x(3) = double(mlon(end,end)); y(3) = double(mlat(end,end));
    x(4) = double(mlon(end,1)); y(4) = double(mlat(end,1));
    x(5) = x(1); y(5) = y(1);
    
    m_line(x,y,'color','r');

    % Put the time in at the middle of the granule
    time_ind = regexp(F(a).name,'\.\d\d\d\d\.');
    time_str = F(a).name(time_ind+1:time_ind+4);
    m_text(mean(x), mean(y), time_str);
end


end

