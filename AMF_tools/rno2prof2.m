%%rno2prof
%%arr 07/24/2008

%..........................................................................
% Returns vector of NO2 mixing ratios on input pressure grid,
% for nearest longitude and latitude in geographic grid
%
% New, simplified version of rCtmGcw (without write capability) that 
% also interpolates profile onto input pressure grid (2008-06-12)
%
% Inputs:
%  fileNO2  = Complete directory/filename of NO2 mixing ratio file
%  pressure = vector of pressures (hPa) monotonically decreasing
%  lon, lat are scalars (deg) as is mon (1,2,...12).
%
% Input keyword:
%  psf = surface pressure (hPa). If present, assume we are reading gmi files.
%        In that case, presSAVE(npres,2) = [[am], [bm]], 
%           where am and bm are vectors or npres elements.
%        Then the model pressure vector is  presModel =  am  +  bm*psf
%
% Outputs (optional: Return these to avoid having to re-read fileNO2 on each call):
%  presSAVE = pressure grid (vector) from fileNO2 (hPA) (unless keyword "psf" is present)
%  lonSAVE = vector of longitudes from fileNO2 (deg)
%  latSAVE = vector of latitudes from fileNO2 (deg)
%  monSAVE = vector of months (1,...12)
%  no2SAVE = large array of no2 mixing ratio values from fileNO2 (npres,nlon,nlat,nmon)
%
%..........................................................................
function [presSAVE, lonSAVE, latSAVE, monSAVE, no2SAVE, no2] = rno2prof2(fileNO2, pressure, lon, lat, mon) %%removed psf bc not reading from gmi files

np = numel(pressure);
dp = circshift(pressure,[1,1])-pressure; 

%if min(dp(1:np-1)) <= 0; 
%    disp('Invalid pressure for rno2prof!')
%else

    %{
    nNO2  = numel(no2SAVE);
    nPres = numel(presSAVE); 
    nLon  = numel(lonSAVE);
    nLat  = numel(latSAVE);  
    nMon  = numel(monSAVE);

    igmi = 0;
    if numel(psf) == 1;
        if psf > 0;
            igmi = 1;
        end
    end
    if igmi == 1;
        nPres = nPres / 2;
    end

    % Read the file (if not already read)
    if nNO2 < 10 || (nNO2 ~= (nPres*nLon*nLat*nMon));
    %}
%fileNO2='C:\Ashley\NewCaRetrieval2\tools\PRFTAV.txt';
if exist('no2SAVE','var')==0;
    
        fid = fopen(fileNO2,'r');

        header = fgets(fid);

        nPres = fscanf(fid,'%g',[1 1]);

        LON = fscanf(fid,'%g',[1 3]);
        nLon = LON(1);  lonStart = LON(2);  lonEnd = LON(3);

        LAT = fscanf(fid,'%g',[1 3]);
        nLat = LAT(1);  latStart = LAT(2);  latEnd = LAT(3);

        MON = fscanf(fid,'%g',[1 3]);
        nMon = MON(1); monStart = MON(2); monEnd = MON(3);
        
        presSAVE = zeros(nPres,1);
        %if igmi == 1;
        %    presSAVE = zeros(nPres,2);  %%rows vs columns issue??
        %end
        presSAVE = fscanf(fid,'%g',[1 nPres]); presSAVE=presSAVE';

        no2SAVE = zeros(nPres,nLon,nLat,nMon);
        no2SAVE_vec = fscanf(fid,'%g',inf);
        no2SAVE = reshape(no2SAVE_vec, [nPres nLon nLat nMon]);
end
        lonSAVE = lonStart + ((0:nLon-1) + 0.5) * (lonEnd - lonStart) / nLon; lonSAVE=lonSAVE';
        latSAVE = latStart + ((0:nLat-1) + 0.5) * (latEnd - latStart) / nLat; latSAVE=latSAVE';
        monSAVE = monStart + ((0:nMon-1) + 0.5) * (monEnd - monStart) / nMon; monSAVE=monSAVE';
    %end
%end

% Interpolate onto input pressure grid

presModel = presSAVE;

for no2_i=1:length(presSAVE)
    dum = reshape(no2SAVE(no2_i,:,:), [nLon, nLat]);
    no21(:,no2_i)=interpn(lonSAVE, latSAVE, dum, lon, lat, 'linear'); 
end
no21(no21==0)=1E-30;
no2=exp(interp1(log(presModel),log(no21'), log(pressure),'linear','extrap'))';

status = fclose(fid);
  