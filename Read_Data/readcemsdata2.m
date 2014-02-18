%readcemsdata2.m
%arr 4/1/2011


addpath C:\Ashley\Trends
addpath J:\CEMS

top200=xlsread('J:\CEMS\top1000_ALLYEARalt.xls','C2:C201');

%years={'2004';'2005';'2006';'2007';'2008';'2009';'2010'};
years={'2011'};
months={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'};
%months={'10';'11';'12'};


for j=1:length(years)
year=years{j};
for k=1:length(months)
month=months{k};  
files=[year,'*.csv'];
cd(['J:/CEMS/',year,'/',month])
cems_files=dir(fullfile('J:','CEMS',year,month,files));
a=0;

howbig=2E7;  % guess how many rows of data and grab the space first

%cemsState=cell(1,howbig); %cemsFacilityName=cell(1,howbig);
cemsSiteCode=zeros(1,howbig); cemsUnitID=cell(1,howbig);
cemsDay=zeros(1,howbig); cemsMonth=zeros(1,howbig);
cemsYear=zeros(1,howbig); cemsHour=zeros(1,howbig);
cemsOpHour=zeros(1,howbig); %cemsGLoad=zeros(1,howbig);
%cemsNOxRate=zeros(1,howbig); cemsNOxRateMeasFlag=cell(1,howbig);
cemsNOxMass=zeros(1,howbig); cemsNOxMassMeasFlag=cell(1,howbig);
    
nlines=1;
for i=1:length(cems_files)
a=a+1;
    if mod(a,20)==0
        disp([num2str(a),' files read so far'])
        %savename=['CEMS',year,'data'];
        %save(savename, 'cemsSiteCode', 'cemsDay', 'cemsMonth', 'cemsYear', 'cemsHour', 'cemsOpHour', 'cemsNOxMass', 'cemsNOxMassMeasFlag');
    end
    fid=fopen(cems_files(i).name, 'r');
    line=fgets(fid);
    line=fgets(fid);
    while line~=-1
    text = textscan(line, '%q', 'delimiter', ',');
    text_char = char(text{1});

    if ismember(str2double(text_char(3,:)),top200)==0;
        line=fgets(fid);
    else
        %{
        thisState=text_char(1,:);
           if isempty(thisState);
              thisState=NaN;
           end
           cemsState{nlines}=thisState;

        thisFacilityName=text_char(2,:);
           if isempty(thisFacilityName);
              thisFacilityName=NaN;
           end
           cemsFacilityName{nlines}=thisFacilityName;
         %}
        thisSiteCode=str2double(text_char(3,:));
           if isempty(thisSiteCode);
              thisSiteCode=NaN;
           end
           cemsSiteCode(nlines)=thisSiteCode;
        
        thisUnitID=text_char(4,:);
           if isempty(thisUnitID);
              thisUnitID=NaN;
           end
           cemsUnitID{nlines}=thisUnitID;
        
        thisDay=str2double(text_char(5,4:5));
           if isempty(thisDay);
              thisDay=NaN;
           end
           cemsDay(nlines)=thisDay;

        thisMonth=str2double(text_char(5,1:2));
           if isempty(thisMonth);
              thisMonth=NaN;
           end
           cemsMonth(nlines)=thisMonth;

        thisYear=str2double(text_char(5,7:10));
           if isempty(thisYear);
              thisYear=NaN;
           end
           cemsYear(nlines)=thisYear;

        thisHour=str2double(text_char(6,:));
           if isempty(thisHour);
              thisHour=NaN;
           end
           cemsHour(nlines)=thisHour;

        thisOpHour=str2double(text_char(7,:));
           if isempty(thisOpHour);
              thisOpHour=NaN;
           end
           cemsOpHour(nlines)=thisOpHour;
        %{
        thisGLoad=str2num(text_char(8,:));
           if isempty(thisGLoad);
              thisGLoad=NaN;
           end
           cemsGLoad(nlines)=thisGLoad;

        thisNOxRate=str2num(text_char(14,:));
           if isempty(thisNOxRate);
              thisNOxRate=NaN;
           end
           cemsNOxRate(nlines)=thisNOxRate;

        thisNOxRateMeasFlag=text_char(15,:);
           if isempty(thisNOxRateMeasFlag);
              thisNOxRateMeasFlag=NaN;
           end
           cemsNOxRateMeasFlag{nlines}=thisNOxRateMeasFlag;
        %}
        thisNOxMass=str2double(text_char(16,:));
           if isempty(thisNOxMass);
              thisNOxMass=NaN;
           end
           cemsNOxMass(nlines)=thisNOxMass;

        thisNOxMassMeasFlag=text_char(17,:);
           if isempty(thisNOxMassMeasFlag);
              thisNOxMassMeasFlag=NaN;
           end
           cemsNOxMassMeasFlag{nlines}=thisNOxMassMeasFlag;
    
        line=fgets(fid);
        nlines=nlines+1;
       %if mod(nlines,2000)==0
       %   disp([num2str(nlines),'lines read so far'])
       %end
    end
    end
fclose(fid); 
end
%cemsState=cemsState(1:nlines); cemsState=deblank(cemsState);
%cemsFacilityName=cemsFacilityName(1:nlines); cemsFacilityName=deblank(cemsFacilityName);
cemsSiteCode=cemsSiteCode(1:nlines-1);
cemsUnitID=cemsUnitID(1:nlines-1); cemsUnitID=deblank(cemsUnitID);
cemsDay=cemsDay(1:nlines-1); %cemsDay=deblank(cemsDay);
cemsMonth=cemsMonth(1:nlines-1); %cemsMonth=deblank(cemsMonth);
cemsYear=cemsYear(1:nlines-1); %cemsYear=deblank(cemsYear);
cemsHour=cemsHour(1:nlines-1);
cemsOpHour=cemsOpHour(1:nlines-1);
%cemsGLoad=cemsGLoad(1:nlines);
%cemsNOxRate=cemsNOxRate(1:nlines);
%cemsNOxRateMeasFlag=cemsNOxRateMeasFlag(1:nlines); cemsNOxRateMeasFlag=deblank(cemsNOxRateMeasFlag);
cemsNOxMass=cemsNOxMass(1:nlines-1);
cemsNOxMassMeasFlag=cemsNOxMassMeasFlag(1:nlines-1); cemsNOxMassMeasFlag=deblank(cemsNOxMassMeasFlag);

       
%{       
cemsWeekday=zeros(1,length(nlines));
for i=1:length(cemsDay);
    day=cemsDay{i};
    month=cemsMonth{i};
    date=[month,'/',day,'/',year];
    a=datenum(date);
    [n, s] = weekday(a);
    cemsWeekday(i)=n;
end
%}
cd(['J:/CEMS/',year])
savename=['CEMS',year,month,'data'];

save(savename, 'cemsUnitID', 'cemsSiteCode', 'cemsDay', 'cemsMonth', 'cemsYear', 'cemsHour', 'cemsOpHour', 'cemsNOxMass', 'cemsNOxMassMeasFlag');
end
end
