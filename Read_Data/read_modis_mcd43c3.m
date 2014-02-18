%read_modis_mcd43c3
%arr 06/21/2010

addpath J:\MODIS_8day

years={'2012'};
months={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'};

Albedo_BSA_Band3=zeros(3600,7200);
Albedo_WSA_Band3=zeros(3600,7200);
BRDF_Quality=zeros(3600,7200);

for k=1:length(years);
    year=years{k};
    file=(['MCD43C3.A',year,'*.hdf']); 
    files=dir(fullfile('J:','MODIS_8day',file));
%count=0;
for i=1:length(months);
    month=months{i};
    count=0;
    for j=1:length(files);
        count=count+1;
        filename=files(count).name;
        hinfo=hdfinfo(filename);  
        
        Albedo_BSA_Band3=hdfread(hinfo.Vgroup(1).Vgroup(1).SDS(3)); Albedo_BSA_Band3=double(Albedo_BSA_Band3);
        Albedo_WSA_Band3=hdfread(hinfo.Vgroup(1).Vgroup(1).SDS(13)); Albedo_WSA_Band3=double(Albedo_WSA_Band3);
        BRDF_Quality=hdfread(hinfo.Vgroup(1).Vgroup(1).SDS(21)); BRDF_Quality=double(BRDF_Quality);

        Albedo_BSA_Band3(Albedo_BSA_Band3==32767)=NaN;
        Albedo_WSA_Band3(Albedo_WSA_Band3==32767)=NaN;
        BRDF_Quality(BRDF_Quality==255)=NaN;
        
        Albedo_BSA_Band3=flipud(Albedo_BSA_Band3); Albedo_WSA_Band3=flipud(Albedo_WSA_Band3); BRDF_Quality=flipud(BRDF_Quality);
        Albedo_BSA_Band3=Albedo_BSA_Band3*0.001; Albedo_WSA_Band3=Albedo_WSA_Band3*0.001;
      
        
        cd J:\MODIS_8day
        n=str2num(filename(14:16));
        date = datestr(addtodate(datenum(['01/01/',year]),n-1, 'day'),2);
        y=['20',date(7:8)];
        m=date(1:2);
        d=date(4:5);
        savename=['mcd43c3_',y,m,d];
        save(savename,'Albedo_BSA_Band3','Albedo_WSA_Band3','BRDF_Quality')
    end
end
end