function [  ] = globe_fix_16aug2016( start_date, end_date )
%GLOBE_FIX_16AUG2016 Fixes an issue with sea-level GLOBE terrain pressure
%   Pixels over ocean in BEHR have a GLOBETerpres value of ~1080 hPa, which
%   is wrong. This happens because a fill value of -500 snuck into the
%   terrain altitude at the beginning wherever GLOBE is a NaN. This code
%   will be just redoing the GLOBE calculation.

DEBUG_LEVEL = 1;

globe_dir = '/global/scratch/laughner/SAT/BEHR/GLOBE_Database/';
load_dir = '/global/scratch/laughner/SAT/BEHR/SP_Files_preGlobeFix';
save_dir = '/global/scratch/laughner/SAT/BEHR/SP_Files';

% These will be set from the data structure
lonlim = [-135 -55];
latlim = [15 60];

% Load the terrain pressure with the new boundaries
[terpres, refvec] = globedem(globe_dir,1,latlim,lonlim);
%refvec will contain (1) number of cells per degree, (2)
%northwest corner latitude, (3) NW corner longitude.
%(2) & (3) might differ from the input latmin & lonmin
%because of where the globe cell edges fall
if DEBUG_LEVEL > 0; fprintf('\n Creating lon/lat matrices for GLOBE data \n'); end
cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1);

terpres(isnan(terpres)) = 0; % CHANGED: NaNs occur over ocean so should be altitude ASL of 0

dnums = datenum(start_date):datenum(end_date);
parfor d=1:numel(dnums)
    tID = getCurrentTask;
    load_name = sprintf('OMI_SP_v2-1Arev1_%04d%02d%02d.mat', year(dnums(d)), month(dnums(d)), day(dnums(d)));
    file_name = fullfile(load_dir, load_name);
    if ~exist(file_name, 'file')
        fprintf('%s does not exist, skipping\n', load_name)
    else
        D = load(file_name);
        if DEBUG_LEVEL > 0; fprintf('w%d: Loaded %s\n', tID.ID, load_name); end
        Data = D.Data;
        for E=1:numel(Data)
            if numel(Data(E).Loncorn) == 1
                fprintf('      w%d: Swath %d has no data, saving as is\n', tID.ID, E)
            else
                if DEBUG_LEVEL > 0; fprintf('    w%d: Swath %d of %d\n', tID.ID, E, numel(Data)); end
                GLOBETerpres = zeros(size(Data(E).Latitude));
                
                %GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
                %corner and terpres(end, end) is in the NE corner.
                
                for k=1:numel(Data(E).Longitude)
                    
                    if DEBUG_LEVEL > 1; fprintf('Averaging GLOBE data to pixel %u of %u \n',k,c); end
                    if DEBUG_LEVEL > 2; tic; end
                    x1 = Data(E).Loncorn(1,k);   y1 = Data(E).Latcorn(1,k);
                    x2 = Data(E).Loncorn(2,k);   y2 = Data(E).Latcorn(2,k);
                    x3 = Data(E).Loncorn(3,k);   y3 = Data(E).Latcorn(3,k);
                    x4 = Data(E).Loncorn(4,k);   y4 = Data(E).Latcorn(4,k);
                    
                    
                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];
                    
                    %%%%SPEED IT UP%%%%
                    % Since GLOBE data is on a grid where a row of
                    % latitudinal points all have the same longitude and
                    % vice versa, we can quickly reduce the number of
                    % points by comparing just one lat and lon vector to
                    % the extent of the pixel.
                    ai=find(globe_lat_matrix(:,1)>=min(yall) & globe_lat_matrix(:,1)<=max(yall));
                    bi=find(globe_lon_matrix(1,:)>=min(xall) & globe_lon_matrix(1,:)<=max(xall));
                    pressurex=terpres(ai,bi);
                    pressure_latx=globe_lat_matrix(ai,bi);
                    pressure_lonx=globe_lon_matrix(ai,bi);
                    %%%%%%%%%%%%%%%%%%%
                    
                    % inpolygon is slow compared to a simple logical test,
                    % so we only apply it to the subset of GLOBE heights
                    % immediately around our pixel.
                    xx_globe = inpolygon(pressure_latx,pressure_lonx,yall,xall);
                    
                    pres_vals=pressurex(xx_globe);
                    GLOBETerpres(k)=1013.25 .* exp(-mean(pres_vals) / 7400 ); %Originally divided by 7640 m
                    
                    
                    if DEBUG_LEVEL > 2; telap = toc; fprintf('Time for GLOBE --> pixel %u/%u = %g sec \n',k,c,telap); end
                end
                
                Data(E).GLOBETerpres = GLOBETerpres;
            end
        end
        save_name = sprintf('OMI_SP_v2-1B_%04d%02d%02d.mat',year(dnums(d)), month(dnums(d)), day(dnums(d)));
        saveData(fullfile(save_dir, save_name),Data);
    end
end


end

function saveData(save_path, Data) %#ok<INUSD>
save(save_path,'Data');
end
