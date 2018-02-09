
newdata = '/Users/monicazhu/Documents/MATLAB/BEHR/BEHR-core/Production tests/UnitTestData/BEHR TropsPres/OMI_BEHR-MONTHLY_US_v3-0A_20120703';
newdir  = '/Users/monicazhu/Documents/MATLAB/BEHR/BEHR-core/Production tests/UnitTestData/BEHR TropsPres';
newpattern = 'OMI_BEHR-MONTHLY_US*';

% generate data;
% read_main('start','2012/07/03','end','2012/07/03','sp_mat_dir',newdir);
% BEHR_main('start','2012/07/03','end','2012/07/03','sp_mat_dir',newdir,'behr_mat_dir',newdir);

new     = load(newdata);
olddata = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v3-0A/DailyProfs/OMI_BEHR_v2-1C_20120703.mat';
olddir  = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v3-0A/DailyProfs';
old     = load(olddata);
%[test]    = behr_unit_test(old.Data,new.Data);
test    = behr_prod_test( newdir,newpattern, olddir  );
%test    = behr_prod_indiv_scatter( newdir, 'OMI_BEHR_v2-1C_20120603.mat', olddir,'OMI_BEHR_v2-1C_20120603.mat'  )

% plot tropopause pressue;
lon = new.OMI.Longitude;
lat = new.OMI.Latitude;
TropoPres = new.OMI.TropopausePressure;

figure;
pcolor(lon,lat,TropoPres);
state_outlines;
colormap(jet);
shading flat;
h=colorbar;