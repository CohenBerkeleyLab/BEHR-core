% readaqsdatainfo.m
%arr 08/09/2011
% data in a string:AA char 1-2, Action code char 3, State Code char 4-5, County Code char 6-8,...
% Site ID char 9-12, Latitude 13-22, Longitude 23-32
cd G:\AQS
tic
filename='Amp500_Site_Monitor-0.txt';
howbig=1.6E7;  % guess how many rows of data and grab the space first
aqsRD=zeros(1,howbig); aqsStateCode=zeros(1,howbig); 
aqsCountyCode=zeros(1,howbig); aqsSiteID=zeros(1,howbig); aqsLatitude=zeros(1,howbig); 
aqsLongitude=zeros(1,howbig); 

fid=fopen(filename, 'r');
line=fgets(fid); line=fgets(fid); line=fgets(fid); line=fgets(fid); 
line=fgets(fid); line=fgets(fid); line=fgets(fid); line=fgets(fid); 
line=fgets(fid); line=fgets(fid); line=fgets(fid); line=fgets(fid);
line=fgets(fid); line=fgets(fid); line=fgets(fid); 
nlines=1;
while 1
   if ~ischar(line), break, end 
   
   text = textscan(line, '%s', 'delimiter', '|');
   text_char = char(text{1});
   %
   thisRD=text_char(1,:); thisRD=deblank(thisRD);
   if isequal(thisRD,'AA')==0
       line=fgets(fid);
   else
   
       thisStateCode=str2double(text_char(3,:));
       if isempty(thisStateCode);
          thisStateCode=NaN;
       end
       aqsStateCode(nlines)=thisStateCode;

       thisCountyCode=str2double(text_char(4,:));
       if isempty(thisCountyCode);
          thisCountyCode=NaN;
       end
       aqsCountyCode(nlines)=thisCountyCode;

       thisSiteID=str2double(text_char(5,:));
       if isempty(thisSiteID);
          thisSiteID=NaN;
       end
       aqsSiteID(nlines)=thisSiteID;

       thisLatitude=str2double(text_char(6,:));
       if isempty(thisLatitude);
          thisLatitude=NaN;
       end
       aqsLatitude(nlines)=thisLatitude;

       thisLongitude=str2double(text_char(7,:));
       if isempty(thisLongitude);
          thisLongitude=NaN;
       end
       aqsLongitude(nlines)=thisLongitude;


       line=fgets(fid);
       nlines=nlines+1;
       if mod(nlines,2000)==0
          disp([num2str(nlines),'lines read so far'])
       end
   end
end
fclose(fid); 


  StateCodes = aqsStateCode(1:nlines);
  CountyCodes = aqsCountyCode(1:nlines);
  SiteIDs = aqsSiteID(1:nlines);
  Latitudes = aqsLatitude(1:nlines);
  Longitudes = aqsLongitude(1:nlines);
  
  
  save AQSdata StateCodes CountyCodes SiteIDs Latitudes Longitudes

toc  