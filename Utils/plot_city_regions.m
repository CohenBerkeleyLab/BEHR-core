%plot_city_regions
%arr 08/12/2011

warning off all

%%%%%%%%%%%%%%SWUS%%%%%%%%%%%%%%%%%%%
%sc_region (south coast)
x=-118.6:0.01:-117.2;
y=sqrt((0.7)^2-(x+117.9).^2)+34;
y2=-sqrt((0.7)^2-(x+117.9).^2)+34;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%fr_region (fresno)
x=-120:0.01:-119.5;
y=sqrt((0.25)^2-(x+119.75).^2)+36.7;
y2=-sqrt((0.25)^2-(x+119.75).^2)+36.7;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%bk_region (bakersfield)
x=-119.25:0.01:-118.75;
y=sqrt((0.25)^2-(x+119).^2)+35.3;
y2=-sqrt((0.25)^2-(x+119).^2)+35.3;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%sd_region (san diego)
x=-117.25:0.01:-116.75;
y=sqrt((0.25)^2-(x+117).^2)+32.8;
y2=-sqrt((0.25)^2-(x+117).^2)+32.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%lv_region (las vegas)
x=-115.45:0.01:-114.95;
y=sqrt((0.25)^2-(x+115.2).^2)+36.2;
y2=-sqrt((0.25)^2-(x+115.2).^2)+36.2;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ph_region (pheonix)
x=-112.4:0.01:-111.6;
y=sqrt((0.4)^2-(x+112).^2)+33.6;
y2=-sqrt((0.4)^2-(x+112).^2)+33.6;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%azpp_region (Power plant in N AZ)
x=-111.55:0.005:-111.15;
y=sqrt((0.2)^2-(x+111.35).^2)+36.9;
y2=-sqrt((0.2)^2-(x+111.35).^2)+36.9;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%fcpp_region (Four Corners power plants)
x=-108.65:0.01:-107.95;
y=sqrt((0.35)^2-(x+108.3).^2)+36.75;
y2=-sqrt((0.35)^2-(x+108.3).^2)+36.75;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%da_region (Dallas)
x=-97.35:0.01:-96.55;
y=sqrt((0.4)^2-(x+96.95).^2)+32.85;
y2=-sqrt((0.4)^2-(x+96.95).^2)+32.85;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%sn_region (San Antonio)
x=-98.75:0.01:-98.15;
y=sqrt((0.3)^2-(x+98.45).^2)+29.55;
y2=-sqrt((0.3)^2-(x+98.45).^2)+29.55;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%al_region (Albequerque)
x=-106.85:0.01:-106.25;
y=sqrt((0.3)^2-(x+106.55).^2)+35.2;
y2=-sqrt((0.3)^2-(x+106.55).^2)+35.2;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%vapp_region (Valmy Power Plant)
x=-117.2:0.01:-116.8;
y=sqrt((0.2)^2-(x+117).^2)+40.9;
y2=-sqrt((0.2)^2-(x+117).^2)+40.9;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%azpp2_region (Arizona Power Plant West)
x=-110.4:0.01:-110;
y=sqrt((0.2)^2-(x+110.2).^2)+34.95;
y2=-sqrt((0.2)^2-(x+110.2).^2)+34.95;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%azpp3_region (Arizona Power Plant East)
x=-109.4:0.01:-109;
y=sqrt((0.2)^2-(x+109.2).^2)+34.55;
y2=-sqrt((0.2)^2-(x+109.2).^2)+34.55;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%tu_region (Tucson)
x=-111.05:0.01:-110.65;
y=sqrt((0.2)^2-(x+110.85).^2)+32.25;
y2=-sqrt((0.2)^2-(x+110.85).^2)+32.25;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%%%%%%%%%%%%%%%NEUS%%%%%%%%%%%%%%%%

%mn_region (minneapolis)
x=-93.55:0.01:-92.95;
y=sqrt((0.3)^2-(x+93.25).^2)+44.95;
y2=-sqrt((0.3)^2-(x+93.25).^2)+44.95;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%kc_region (kansas city)
x=-94.85:0.01:-94.25;
y=sqrt((0.3)^2-(x+94.55).^2)+39.15;
y2=-sqrt((0.3)^2-(x+94.55).^2)+39.15;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%sl_region (st louis)
x=-90.65:0.01:-90.05;
y=sqrt((0.3)^2-(x+90.35).^2)+38.65;
y2=-sqrt((0.3)^2-(x+90.35).^2)+38.65;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ch_region (chicago)
x=-88.3:0.01:-87.1;
y=sqrt((0.6)^2-(x+87.7).^2)+41.8;
y2=-sqrt((0.6)^2-(x+87.7).^2)+41.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%in_region (indianapolis)
x=-86.35:0.01:-85.95;
y=sqrt((0.2)^2-(x+86.15).^2)+39.8;
y2=-sqrt((0.2)^2-(x+86.15).^2)+39.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%de_region (detroit)
x=-83.55:0.01:-82.65;
y=sqrt((0.45)^2-(x+83.1).^2)+42.35;
y2=-sqrt((0.45)^2-(x+83.1).^2)+42.35;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%cl_region (cleveland)
x=-81.97:0.01:-81.37;
y=sqrt((0.3)^2-(x+81.67).^2)+41.45;
y2=-sqrt((0.3)^2-(x+81.67).^2)+41.45;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%co_region (columbus)
x=-83.4:0.01:-82.8;
y=sqrt((0.3)^2-(x+83.1).^2)+40;
y2=-sqrt((0.3)^2-(x+83.1).^2)+40;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%to_region (toronto)
x=-79.9:0.01:-79.1;
y=sqrt((0.4)^2-(x+79.5).^2)+43.7;
y2=-sqrt((0.4)^2-(x+79.5).^2)+43.7;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%mo_region (montreal)
x=-74:0.01:-73.4;
y=sqrt((0.3)^2-(x+73.7).^2)+45.6;
y2=-sqrt((0.3)^2-(x+73.7).^2)+45.6;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%pi_region (pittsburg)
x=-80.2:.01:-79.7;
y=sqrt((0.25)^2-(x+79.95).^2)+40.4;
y2=-sqrt((0.25)^2-(x+79.95).^2)+40.4;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%bo_region (boston)
x=-71.3:.01:-70.7;
y=sqrt((0.3)^2-(x+71).^2)+42.45;
y2=-sqrt((0.3)^2-(x+71).^2)+42.45;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%dc_region (washington dc)
x=-77.3:.01:-76.3;
y=sqrt((0.5)^2-(x+76.8).^2)+39.15;
y2=-sqrt((0.5)^2-(x+76.8).^2)+39.15;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ny_region (new york)
x=-74.4:.01:-73;
y=sqrt((0.7)^2-(x+73.7).^2)+40.85;
y2=-sqrt((0.7)^2-(x+73.7).^2)+40.85;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ph_region (philidelphia)
x=-75.7:.01:-74.7;
y=sqrt((0.5)^2-(x+75.2).^2)+40;
y2=-sqrt((0.5)^2-(x+75.2).^2)+40;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ci_region (cincinnati)
x=-84.8:.01:-84.3;
y=sqrt((0.25)^2-(x+84.55).^2)+39.1;
y2=-sqrt((0.25)^2-(x+84.55).^2)+39.1;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%inpp_region (Indiana Power Plant)
x=-87.2:.01:-86.8;
y=sqrt((0.2)^2-(x+87).^2)+37.95;
y2=-sqrt((0.2)^2-(x+87).^2)+37.95;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ohpp_region (Ohio Power Plant North)
x=-80.65:.01:-80.25;
y=sqrt((0.2)^2-(x+80.45).^2)+40.55;
y2=-sqrt((0.2)^2-(x+80.45).^2)+40.55;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ohpp2_region (Ohio Power Plant South)
x=-80.9:.01:-80.5;
y=sqrt((0.2)^2-(x+80.7).^2)+39.85;
y2=-sqrt((0.2)^2-(x+80.7).^2)+39.85;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%papp_region (Pennsylvania Plant North West)
x=-79.25:.01:-78.85;
y=sqrt((0.2)^2-(x+79.05).^2)+40.4;
y2=-sqrt((0.2)^2-(x+79.05).^2)+40.4;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%papp2_region (Pennsylvania Power Plant South West)
x=-80.1:.01:-79.7;
y=sqrt((0.2)^2-(x+79.9).^2)+39.85;
y2=-sqrt((0.2)^2-(x+79.9).^2)+39.85;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%papp3_region (Pennsylvania Power Plant East)
x=-76.75:.01:-76.35;
y=sqrt((0.2)^2-(x+76.55).^2)+40.1;
y2=-sqrt((0.2)^2-(x+76.55).^2)+40.1;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ri_region (Richmond, VA)
x=-77.5:.01:-77.1;
y=sqrt((0.2)^2-(x+77.3).^2)+37.4;
y2=-sqrt((0.2)^2-(x+77.3).^2)+37.4;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%%%%%%%%%%%%%NWUS%%%%%%%%%%%%%%%%%

%sa_region (sacramento)
x=-121.65:.01:-121.15;
y=sqrt((0.25)^2-(x+121.4).^2)+38.65;
y2=-sqrt((0.25)^2-(x+121.4).^2)+38.65;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%po_region (portland)
x=-122.85:.01:-122.25;
y=sqrt((0.3)^2-(x+122.55).^2)+45.45;
y2=-sqrt((0.3)^2-(x+122.55).^2)+45.45;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%se_region (seattle)
x=-122.75:.01:-121.75;
y=sqrt((0.50)^2-(x+122.25).^2)+47.35;
y2=-sqrt((0.50)^2-(x+122.25).^2)+47.35;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%va_region (vancouver)
x=-123.1:.01:-122.6;
y=sqrt((0.25)^2-(x+122.85).^2)+49.25;
y2=-sqrt((0.25)^2-(x+122.85).^2)+49.25;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%sl_region (salt lake city)
x=-112.15:.01:-111.75;
y=sqrt((0.2)^2-(x+111.95).^2)+40.7;
y2=-sqrt((0.2)^2-(x+111.95).^2)+40.7;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%utpp1_region (utah power plant left)
x=-112.75:.01:-112.35;
y=sqrt((0.2)^2-(x+112.55).^2)+39.55;
y2=-sqrt((0.2)^2-(x+112.55).^2)+39.55;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%utpp2_region (utah power plant right)
x=-111.3:.01:-110.9;
y=sqrt((0.2)^2-(x+111.1).^2)+39.25;
y2=-sqrt((0.2)^2-(x+111.1).^2)+39.25;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%mtpp_region (montana power plant)
x=-106.65:.01:-106.25;
y=sqrt((0.2)^2-(x+106.45).^2)+45.9;
y2=-sqrt((0.2)^2-(x+106.45).^2)+45.9;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%re_region (reno)
x=-119.9:.01:-119.5;
y=sqrt((0.2)^2-(x+119.7).^2)+39.55;
y2=-sqrt((0.2)^2-(x+119.7).^2)+39.55;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%wypp1_region (wyoming power plant, left)
x=-108.85:.01:-108.45;
y=sqrt((0.2)^2-(x+108.65).^2)+41.75;
y2=-sqrt((0.2)^2-(x+108.65).^2)+41.75;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%wypp2_region (wyoming power plant, right)
x=-105.05:.01:-104.65;
y=sqrt((0.2)^2-(x+104.85).^2)+42.1;
y2=-sqrt((0.2)^2-(x+104.85).^2)+42.1;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%dn_region (denver)
x=-105.3:.01:-104.7;
y=sqrt((0.3)^2-(x+105).^2)+39.75;
y2=-sqrt((0.3)^2-(x+105).^2)+39.75;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%copp_region (colorado power plant)
x=-107.55:.01:-107.15;
y=sqrt((0.2)^2-(x+107.35).^2)+40.5;
y2=-sqrt((0.2)^2-(x+107.35).^2)+40.5;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%om_region (omaha)
x=-96.3:.01:-95.8;
y=sqrt((0.25)^2-(x+96.05).^2)+41.3;
y2=-sqrt((0.25)^2-(x+96.05).^2)+41.3;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ndpp_region (north dakota power plants)
x=-101.7:.01:-101;
y=sqrt((0.35)^2-(x+101.35).^2)+47.4;
y2=-sqrt((0.35)^2-(x+101.35).^2)+47.4;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)


%%%%%%%%%%%%%%%%SEUS%%%%%%%%%%%%%%%%%%%

%at_region (atlanta)
x=-84.7:.01:-84;
y=sqrt((0.35)^2-(x+84.35).^2)+33.8;
y2=-sqrt((0.35)^2-(x+84.35).^2)+33.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%no_region (new orleans)
x=-90.6:.01:-90;
y=sqrt((0.3)^2-(x+90.3).^2)+30.05;
y2=-sqrt((0.3)^2-(x+90.3).^2)+30.05;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%tm_region (tampa)
x=-82.7:.01:-82.1;
y=sqrt((0.3)^2-(x+82.4).^2)+27.9;
y2=-sqrt((0.3)^2-(x+82.4).^2)+27.9;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%jk_region (jacksonville)
x=-81.85:.01:-81.35;
y=sqrt((0.25)^2-(x+81.6).^2)+30.45;
y2=-sqrt((0.25)^2-(x+81.6).^2)+30.45;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ch_region (charlotte)
x=-81.05:.01:-80.65;
y=sqrt((0.2)^2-(x+80.85).^2)+35.25;
y2=-sqrt((0.2)^2-(x+80.85).^2)+35.25;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%flpp_region (florida powerplant)
x=-82.85:.01:-82.45;
y=sqrt((0.2)^2-(x+82.65).^2)+28.95;
y2=-sqrt((0.2)^2-(x+82.65).^2)+28.95;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%me_region (memphis)
x=-90.3:.01:-89.9;
y=sqrt((0.2)^2-(x+90.1).^2)+35.1;
y2=-sqrt((0.2)^2-(x+90.1).^2)+35.1;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%na_region (nashville)
x=-86.9:.01:-86.3;
y=sqrt((0.3)^2-(x+86.6).^2)+36.2;
y2=-sqrt((0.3)^2-(x+86.6).^2)+36.2;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%tnpp_region (tennessee power plant)
x=-88.1:.01:-87.7;
y=sqrt((0.2)^2-(x+87.9).^2)+36.05;
y2=-sqrt((0.2)^2-(x+87.9).^2)+36.05;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%kn_region (knoxville)
x=-84.3:.01:-83.7;
y=sqrt((0.3)^2-(x+84).^2)+35.95;
y2=-sqrt((0.3)^2-(x+84).^2)+35.95;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%or_region (orlando)
x=-81.5:.01:-81.1;
y=sqrt((0.2)^2-(x+81.3).^2)+28.5;
y2=-sqrt((0.2)^2-(x+81.3).^2)+28.5;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%mi_region (miami)
x=-80.5:.01:-80.1;
y=sqrt((0.2)^2-(x+80.3).^2)+26.05;
y2=-sqrt((0.2)^2-(x+80.3).^2)+26.05;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%flpp2_region (Florida Power Plant East)
x=-81.75:.01:-81.35;
y=sqrt((0.2)^2-(x+81.55).^2)+29.8;
y2=-sqrt((0.2)^2-(x+81.55).^2)+29.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ncpp_region (North Carolina Power Plant)
x=-81.15:.01:-80.75;
y=sqrt((0.2)^2-(x+80.95).^2)+35.6;
y2=-sqrt((0.2)^2-(x+80.95).^2)+35.6;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ilpp_region (Illinois Power Plant)
x=-88.95:.01:-88.55;
y=sqrt((0.2)^2-(x+88.75).^2)+37.2;
y2=-sqrt((0.2)^2-(x+88.75).^2)+37.2;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%%%%%%%%%%%%%%%%%US%%%%%%%%%%%%%%%%%%

%sf_region (San Francisco)
x=-122.4:.01:-121.6;
y=sqrt((0.4)^2-(x+122).^2)+37.6;
y2=-sqrt((0.4)^2-(x+122).^2)+37.6;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)

%ho_region (Houston)
x=-95.55:.01:-94.95;
y=sqrt((0.3)^2-(x+95.25).^2)+29.8;
y2=-sqrt((0.3)^2-(x+95.25).^2)+29.8;
hold on
plot(x,y,'k-','LineWidth',1)
hold on
plot(x,y2,'k-','LineWidth',1)


%%%%%%%%%%%%%%%%Background%%%%%%%%%%%%%%%%%%%%%

%mt_bkgd (Background taken from Western Montana)
hold on
x=-113:0.1:-108; y=45.5:0.1:48.5; 
plot(x,ones(1,numel(x))*45.5,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*48.5,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-113,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-108,y,'k-','LineWidth',1)


%mt_bkgd2 (Background taken from Northeastern Montana)
hold on
x=-108:0.1:-104.5; y=46.5:0.1:48.5; 
plot(x,ones(1,numel(x))*46.5,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*48.5,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-108,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-104.5,y,'k-','LineWidth',1)

%wy_bkgd (Background taken from Western Wyoming)
hold on
x=-111:0.1:-107; y=42.5:0.1:45; 
plot(x,ones(1,numel(x))*42.5,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*45,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-111,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-107,y,'k-','LineWidth',1)

%or_bkgd (Background taken from Southeastern Oregon)
hold on
x=-121:0.1:-117; y=42:0.1:44.5; 
plot(x,ones(1,numel(x))*42,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*44.5,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-121,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-117,y,'k-','LineWidth',1)

%nv_bkgd (Background taken from central Nevada)
hold on
x=-118:0.1:-114.5; y=38:0.1:40; 
plot(x,ones(1,numel(x))*38,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*40,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-118,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-114.5,y,'k-','LineWidth',1)

%wa_bkgd (Background taken from eastern washington)
hold on
x=-121:0.1:-117; y=48:0.1:49; 
plot(x,ones(1,numel(x))*48,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*49,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-121,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-117,y,'k-','LineWidth',1)

%sd_bkgd (Background taken from south dakota)
hold on
x=-104:0.1:-100; y=43:0.1:46; 
plot(x,ones(1,numel(x))*43,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*46,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-104,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-100,y,'k-','LineWidth',1)


%me_bkgd (Background taken from Maine)
hold on
x=-70:0.1:-68; y=45:0.1:46.5; 
plot(x,ones(1,numel(x))*45,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*46.5,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-70,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-68,y,'k-','LineWidth',1)

%in_bkgd (Background taken from northern indiana)
hold on
x=-88.5:0.1:-84.5; y=40.2:0.1:41; 
plot(x,ones(1,numel(x))*40.2,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*41,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-88.5,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-84.5,y,'k-','LineWidth',1)

%tx_bkgd (Background taken from western Texas)
hold on
x=-102.5:0.1:-99; y=30:0.1:33; 
plot(x,ones(1,numel(x))*30,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*33,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-102.5,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-99,y,'k-','LineWidth',1)

%nm_bkgd (Background taken from New Mexico)
hold on
x=-108.5:0.1:-105; y=32.5:0.1:34.5; 
plot(x,ones(1,numel(x))*32.5,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*34.5,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-108.5,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-105,y,'k-','LineWidth',1)

%pac_bkgd (Background taken from the Pacific off SoCal)
hold on
x=-124:0.1:-121; y=32:0.1:35; 
plot(x,ones(1,numel(x))*32,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*35,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-124,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-121,y,'k-','LineWidth',1)

%mex_bkgd (Background taken from northern mexico)
hold on
x=-111:0.1:-104; y=30:0.1:31; 
plot(x,ones(1,numel(x))*30,'k-','LineWidth',1)
hold on
plot(x,ones(1,numel(x))*31,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-111,y,'k-','LineWidth',1)
hold on
plot(ones(1,numel(y))*-104,y,'k-','LineWidth',1)

