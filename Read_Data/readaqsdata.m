% readaqsdata.m
% read 2007 aqs data from epa
%arr 05/18/2010
% data in a string:RD char 1-2, Action code char 3, State Code char 4-5, County Code char 6-8,...
% Site ID char 9-12, Parameter char 13-17, POC char 18, Sample duration char 19, Unit char 20-22,...
% Method char 23-25, Date char 26-33, Start Time char 34-38, Sample Value char 64;

cd G:\AQS
tic
filename='RD_501_42602_2010-0.txt';
howbig=1.6E7;  % guess how many rows of data and grab the space first
aqsRD=zeros(1,howbig); aqsActionCode=zeros(1,howbig); aqsStateCode=zeros(1,howbig); 
aqsCountyCode=zeros(1,howbig); aqsSiteID=zeros(1,howbig); aqsParameter=zeros(1,howbig); 
aqsPOC=zeros(1,howbig); aqsSampleDuration=zeros(1,howbig); aqsUnit=zeros(1,howbig); 
aqsMethod=zeros(1,howbig); aqsYear=zeros(1,howbig); aqsMonth=zeros(1,howbig);
aqsDay=zeros(1,howbig); aqsHour=zeros(1,howbig); aqsNO2=zeros(1,howbig);

fid=fopen(filename, 'r');
line=fgets(fid);
line=fgets(fid);
line=fgets(fid);
nlines=1;
while 1
   if ~isstr(line), break, end 
   
   text = textscan(line, '%s', 'delimiter', '|');
   text_char = char(text{1});
   %{
   thisRD=text_char(1,:);
   if isempty(thisRD);
       thisRD=NaN;
   end
   aqsRD(nlines)=thisRD;
   
   thisActionCode=text_char(2,:);
   if isempty(thisActionCode);
      thisActionCode=NaN;
   end
   aqsActionCode(nlines)=thisActionCode;
   %}
   
   thisStateCode=str2num(text_char(3,:));
   if isempty(thisStateCode);
      thisStateCode=NaN;
   end
   aqsStateCode(nlines)=thisStateCode;
   
   thisCountyCode=str2num(text_char(4,:));
   if isempty(thisCountyCode);
      thisCountyCode=NaN;
   end
   aqsCountyCode(nlines)=thisCountyCode;
   
   thisSiteID=str2num(text_char(5,:));
   if isempty(thisSiteID);
      thisSiteID=NaN;
   end
   aqsSiteID(nlines)=thisSiteID;
   
   thisParameter=str2num(text_char(6,:));
   if isempty(thisParameter);
      thisParameter=NaN;
   end
   aqsParameter(nlines)=thisParameter;
   
   thisPOC=str2num(text_char(7,:));
   if isempty(thisPOC);
      thisPOC=NaN;
   end
   aqsPOC(nlines)=thisPOC;
   
   thisSampleDuration=str2num(text_char(8,:));
   if isempty(thisSampleDuration);
      thisSampleDuration=NaN;
   end
   aqsSampleDuration(nlines)=thisSampleDuration;
   
   thisUnit=str2num(text_char(9,:));
   if isempty(thisUnit);
      thisUnit=NaN;
   end
   aqsUnit(nlines)=thisUnit;
   
   thisMethod=str2num(text_char(10,:));
   if isempty(thisMethod);
      thisMethod=NaN;
   end
   aqsMethod(nlines)=thisMethod;
   
   thisYear=str2num(text_char(11,1:4));
   if isempty(thisYear);
      thisYear=NaN;
   end
   aqsYear(nlines)=thisYear;
   
   thisMonth=str2num(text_char(11,5:6));
   if isempty(thisMonth);
      thisMonth=NaN;
   end
   aqsMonth(nlines)=thisMonth;
   
   thisDay=str2num(text_char(11,7:8));
   if isempty(thisDay);
      thisDay=NaN;
   end
   aqsDay(nlines)=thisDay;
   
   thisHour=str2num(text_char(12,1:2));
   if isempty(thisHour);
      thisHour=NaN;
   end
   aqsHour(nlines)=thisHour;
   
   thisNO2=str2num(text_char(13,:));
   if isempty(thisNO2);
      thisNO2=NaN;
   end
   aqsNO2(nlines)=thisNO2;
   
   
   line=fgets(fid);
   nlines=nlines+1;
   if mod(nlines,2000)==0
      disp([num2str(nlines),'lines read so far'])
   end
   % end
end
fclose(fid); 


  aqsRD = aqsRD(1:nlines);
  aqsActionCode = aqsActionCode(1:nlines);
  aqsStateCode = aqsStateCode(1:nlines);
  aqsCountyCode = aqsCountyCode(1:nlines);
  aqsSiteID = aqsSiteID(1:nlines);
  aqsParameter = aqsParameter(1:nlines);
  aqsPOC = aqsPOC(1:nlines);
  aqsSampleDuration = aqsSampleDuration(1:nlines);
  aqsUnit = aqsUnit(1:nlines);
  aqsMethod =aqsMethod(1:nlines);
  aqsYear = aqsYear(1:nlines);
  aqsMonth = aqsMonth(1:nlines);
  aqsDay = aqsDay(1:nlines);
  aqsHour = aqsHour(1:nlines);
  aqsNO2 = aqsNO2(1:nlines);
  
  save AQS2010data aqsRD aqsActionCode aqsStateCode aqsCountyCode aqsSiteID aqsParameter aqsPOC aqsSampleDuration aqsUnit aqsMethod aqsYear aqsMonth aqsDay aqsHour aqsNO2

toc  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add Latitude and Longitude to AQS(year)data mat-files

years={'2004';'2005';'2006';'2007';'2008';'2009';'2010'};
for i=4;%:length(years);
    year=years{i};
load(['AQS',year,'dataO3'])
load('G:\AQS\AQSdata.mat')

aqsLatitude=zeros(1,length(aqsDay));
aqsLongitude=zeros(1,length(aqsDay));
for i=1:length(aqsDay)
    a=find(StateCodes==aqsStateCode(i) & CountyCodes==aqsCountyCode(i) & SiteIDs==aqsSiteID(i));
    if isempty(a)
    else
        aqsLatitude(i)=Latitudes(a);
        aqsLongitude(i)=Longitudes(a);
    end
end

savename=['AQS',year,'dataO3'];
cd G:\AQS\
save(savename, 'aqsRD', 'aqsActionCode', 'aqsStateCode', 'aqsCountyCode', 'aqsSiteID', 'aqsParameter', 'aqsPOC', 'aqsSampleDuration', 'aqsUnit', 'aqsMethod', 'aqsYear', 'aqsMonth', 'aqsDay', 'aqsHour', 'aqsO3', 'aqsLatitude', 'aqsLongitude');
end