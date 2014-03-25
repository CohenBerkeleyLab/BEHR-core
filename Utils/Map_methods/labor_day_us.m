%%labor_day_us
%%arr 08/16/2011

addpath C:\Ashley
addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map1.4
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map
addpath C:\Ashley\Old_stuff
load outside

tic
satellite='OMI';
retrieval='NEW';

load sf_region_us
xx=sf_region_us;

%year_i={'2005';'2005';'2005';'2005';'2006';'2006';'2006';'2006';'2007';'2007';'2007';'2007';'2008';'2008';'2008';'2008';'2009';'2009';'2009';'2009';'2010';'2010';'2010';'2010';'2011';'2011';'2011';'2011'};
%season_i={'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer';'winter';'winter';'summer';'summer'};
%wdwe_i={'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we'};


%year_i={'2005';'2005';'2006';'2006';'2007';'2007';'2008';'2008';'2009';'2009';'2010';'2010';'2011';'2011'};
%season_i={'winter';'summer';'winter';'summer';'winter';'summer';'winter';'summer';'winter';'summer';'winter';'summer';'winter';'summer'};
%wdwe_i={'all';'all';'all';'all';'all';'all';'all';'all';'all';'all';'all';'all';'all';'all'};

year_i={'2009';'2009';'2010';'2010';'2011';'2011'};
season_i={'winter';'summer';'winter';'summer';'winter';'summer'};
wdwe_i={'wd';'wd';'wd';'wd';'wd';'wd'};

%year_i={'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2011';'2011';'2011';'2011';'2011';'2011'};
%season_i={'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum';'fal';'fal';'win';'win';'spr';'spr';'sum';'sum'};
%wdwe_i={'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we';'wd';'we'};



lats=25:0.05:50;
lons=-125:0.05:-65;

allDATA=[];
ALLDATA=[];
index=0;
for x=1:length(year_i);
    index=index+1;
    year=year_i{x};
    season=season_i{x};
    wdwe=wdwe_i{x};

    omi_no_holiday
    
    SumWeightedDataNWUS=zeros(250,600);
    SumWeightedDataSWUS=zeros(250,600);
    SumWeightedDataNEUS=zeros(250,600);
    SumWeightedDataSEUS=zeros(250,600);
    SumWeightNWUS=zeros(250,600);
    SumWeightSWUS=zeros(250,600);
    SumWeightNEUS=zeros(250,600);
    SumWeightSEUS=zeros(250,600);
    
    L=length(days);
    for i=1:L;
        year=years{i};
        month=months{i};
        day=days{i};
        filename=(['OMI_NEW_',year,month,day]); 
        warning off all
        
        %NWUS
        cd('J:\New_Retrieval_sp_mats\NWUS_wCld\BEHR_new')
        a=exist([filename,'.mat'],'file');
        if a==0;
            disp(['No Data Available For ',month,' ',day,' ',year])
        else
            load(filename)
            s=size(OMI);
            for ss=1:s(2);
                omi=OMI(ss);
                aa=find(omi.NewColumnAmountNO2Polluted<=0);       omi.Areaweight(aa)=0;
                %bb=find(omi.vcdQualityFlags==0);               omi.Areaweight(bb)=0;  %%%flag SP=0; NRT=-1;
                bb=find(omi.vcdQualityFlags/2~=round(omi.vcdQualityFlags/2));      omi.Areaweight(bb)=0;
                cc=find(omi.CloudFraction>200);                omi.Areaweight(cc)=0;
                gg=find(omi.NewColumnAmountNO2Polluted>1E17); omi.Areaweight(gg)=0;%omi.NewColumnAmountNO2Polluted(dd)=0;
                hh=find(isnan(omi.NewColumnAmountNO2Polluted)==1); omi.NewColumnAmountNO2Polluted(hh)=0; omi.Areaweight(hh)=0;
                
                %%only throw out bad pixels; leave edge pixels
                %dd=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>0 & omi.ViewingZenithAngle<30); omi.Areaweight(dd)=0;
                %ee=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>48 & omi.ViewingZenithAngle<52); omi.Areaweight(ee)=0;
                %ff=find(omi.ViewingAzimuthAngle>0 & omi.ViewingZenithAngle<7); omi.Areaweight(ff)=0;
                
                %throw out bad pixels using new row variable
                dd=find(omi.Row>=27 & omi.Row<=45); omi.Areaweight(dd)=0;
                ee=find(omi.Row>=53 & omi.Row<=54); omi.Areaweight(ee)=0;
                
                %dd=find(omi.ViewingZenithAngle>22); omi.Areaweight(dd)=0;
                SumWeightedDataNWUS=SumWeightedDataNWUS+(omi.NewColumnAmountNO2Polluted).*omi.Areaweight;
                SumWeightNWUS=SumWeightNWUS+omi.Areaweight;
            end
        end 
        
        %SWUS
        cd('J:\New_Retrieval_sp_mats\SWUS_wCld\BEHR_new')
        a=exist([filename,'.mat'],'file');
        if a==0;
            disp(['No Data Available For ',month,' ',day,' ',year])
        else
            load(filename)
            s=size(OMI);
            for ss=1:s(2);
                omi=OMI(ss);
                aa=find(omi.NewColumnAmountNO2Polluted<=0);       omi.Areaweight(aa)=0;
                %bb=find(omi.vcdQualityFlags==0);               omi.Areaweight(bb)=0;  %%%flag SP=0; NRT=-1;
                bb=find(omi.vcdQualityFlags/2~=round(omi.vcdQualityFlags/2));      omi.Areaweight(bb)=0;
                cc=find(omi.CloudFraction>200);                omi.Areaweight(cc)=0;
                gg=find(omi.NewColumnAmountNO2Polluted>1E17); omi.Areaweight(gg)=0;%omi.NewColumnAmountNO2Polluted(dd)=0;
                hh=find(isnan(omi.NewColumnAmountNO2Polluted)==1); omi.NewColumnAmountNO2Polluted(hh)=0; omi.Areaweight(hh)=0;
                
                %%only throw out bad pixels; leave edge pixels
                %dd=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>0 & omi.ViewingZenithAngle<30); omi.Areaweight(dd)=0;
                %ee=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>48 & omi.ViewingZenithAngle<52); omi.Areaweight(ee)=0;
                %ff=find(omi.ViewingAzimuthAngle>0 & omi.ViewingZenithAngle<7); omi.Areaweight(ff)=0;
                
                %throw out bad pixels using new row variable
                dd=find(omi.Row>=27 & omi.Row<=45); omi.Areaweight(dd)=0;
                ee=find(omi.Row>=53 & omi.Row<=54); omi.Areaweight(ee)=0;
                
                %dd=find(omi.ViewingZenithAngle>22); omi.Areaweight(dd)=0;
                SumWeightedDataSWUS=SumWeightedDataSWUS+(omi.NewColumnAmountNO2Polluted).*omi.Areaweight;
                SumWeightSWUS=SumWeightSWUS+omi.Areaweight;
            end
        end 
        
        %NEUS
        cd('J:\New_Retrieval_sp_mats\NEUS_wCld\BEHR_new')
        a=exist([filename,'.mat'],'file');
        if a==0;
            disp(['No Data Available For ',month,' ',day,' ',year])
        else
            load(filename)
            s=size(OMI);
            for ss=1:s(2);
                omi=OMI(ss);
                aa=find(omi.NewColumnAmountNO2Polluted<=0);       omi.Areaweight(aa)=0;
                %bb=find(omi.vcdQualityFlags==0);               omi.Areaweight(bb)=0;  %%%flag SP=0; NRT=-1;
                bb=find(omi.vcdQualityFlags/2~=round(omi.vcdQualityFlags/2));      omi.Areaweight(bb)=0;
                cc=find(omi.CloudFraction>200);                omi.Areaweight(cc)=0;
                gg=find(omi.NewColumnAmountNO2Polluted>1E17); omi.Areaweight(gg)=0;%omi.NewColumnAmountNO2Polluted(dd)=0;
                hh=find(isnan(omi.NewColumnAmountNO2Polluted)==1); omi.NewColumnAmountNO2Polluted(hh)=0; omi.Areaweight(hh)=0;
                
                %%only throw out bad pixels; leave edge pixels
                %dd=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>0 & omi.ViewingZenithAngle<30); omi.Areaweight(dd)=0;
                %ee=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>48 & omi.ViewingZenithAngle<52); omi.Areaweight(ee)=0;
                %ff=find(omi.ViewingAzimuthAngle>0 & omi.ViewingZenithAngle<7); omi.Areaweight(ff)=0;
                
                %throw out bad pixels using new row variable
                dd=find(omi.Row>=27 & omi.Row<=45); omi.Areaweight(dd)=0;
                ee=find(omi.Row>=53 & omi.Row<=54); omi.Areaweight(ee)=0;
                
                %dd=find(omi.ViewingZenithAngle>22); omi.Areaweight(dd)=0;
                SumWeightedDataNEUS=SumWeightedDataNEUS+(omi.NewColumnAmountNO2Polluted).*omi.Areaweight;
                SumWeightNEUS=SumWeightNEUS+omi.Areaweight;
            end
        end 
        
        %SEUS
        cd('J:\New_Retrieval_sp_mats\SEUS_wCld\BEHR_new')
        a=exist([filename,'.mat'],'file');
        if a==0;
            disp(['No Data Available For ',month,' ',day,' ',year])
        else
            load(filename)
            s=size(OMI);
            for ss=1:s(2);
                omi=OMI(ss);
                aa=find(omi.NewColumnAmountNO2Polluted<=0);       omi.Areaweight(aa)=0;
                %bb=find(omi.vcdQualityFlags==0);               omi.Areaweight(bb)=0;  %%%flag SP=0; NRT=-1;
                bb=find(omi.vcdQualityFlags/2~=round(omi.vcdQualityFlags/2));      omi.Areaweight(bb)=0;
                cc=find(omi.CloudFraction>200);                omi.Areaweight(cc)=0;
                gg=find(omi.NewColumnAmountNO2Polluted>1E17); omi.Areaweight(gg)=0;%omi.NewColumnAmountNO2Polluted(dd)=0;
                hh=find(isnan(omi.NewColumnAmountNO2Polluted)==1); omi.NewColumnAmountNO2Polluted(hh)=0; omi.Areaweight(hh)=0;
                
                %%only throw out bad pixels; leave edge pixels
                %dd=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>0 & omi.ViewingZenithAngle<30); omi.Areaweight(dd)=0;
                %ee=find(omi.ViewingAzimuthAngle<0 & omi.ViewingZenithAngle>48 & omi.ViewingZenithAngle<52); omi.Areaweight(ee)=0;
                %ff=find(omi.ViewingAzimuthAngle>0 & omi.ViewingZenithAngle<7); omi.Areaweight(ff)=0;
                
                %throw out bad pixels using new row variable
                dd=find(omi.Row>=27 & omi.Row<=45); omi.Areaweight(dd)=0;
                ee=find(omi.Row>=53 & omi.Row<=54); omi.Areaweight(ee)=0;
                
              
                %dd=find(omi.ViewingZenithAngle>22); omi.Areaweight(dd)=0;
                SumWeightedDataSEUS=SumWeightedDataSEUS+(omi.NewColumnAmountNO2Polluted).*omi.Areaweight;
                SumWeightSEUS=SumWeightSEUS+omi.Areaweight;
            end
        end 
        
    end
    DataNWUS=SumWeightedDataNWUS./SumWeightNWUS; 
    DataSWUS=SumWeightedDataSWUS./SumWeightSWUS; 
    DataNEUS=SumWeightedDataNEUS./SumWeightNEUS; 
    DataSEUS=SumWeightedDataSEUS./SumWeightSEUS; 

    Data=[DataSWUS DataSEUS;DataNWUS DataNEUS];
    
  
    a=find(isnan(Data)==1); Data(a)=[]; 
    lat=25.025:0.05:49.975; lat=lat'; lat=repmat(lat,1,1200); latx=lat; latx(a)=[];
    lon=-124.975:0.05:-65.025; lon=repmat(lon,500,1); lonx=lon; lonx(a)=[];

    X=griddata(lonx,latx,Data,lon,lat);
    X(numel(lats),:)=0; X(:,numel(lons))=0;
    X(outside==0)=NaN; 
%
    DATA=X.*xx;
    DATA=DATA(:);
    DATA(isnan(DATA))=[]; DATA(DATA==0)=[];
    DATA=mean(DATA)
    allDATA=[allDATA;DATA];
    
    ALLDATA(index,:,:)=X;
%
    
end

%
for i=1:251,
    for j=1:601;
        a=ALLDATA(:,i,j);
        a(isnan(a))=[];
        b=mean(a);
        alldata(i,j)=b;
    end
end
%

figure, m_proj('Mercator', 'long', [-125 -65], 'lat', [25 50]);
m_pcolor(lons,lats,X); shading('flat')
DD=[0 1E16]; CAXIS(DD); colorbar; hold on
plot_states,
m_grid('linestyle','none');
set(gca,'DataAspectRatio',[1 0.8 1])