function [  ] = no2_daily_column_map( date_in, lonlim, latlim, varargin )
%NO2_DAILY_COLUMN_MAP Plots all swaths on a given day inside geo limits

E = JLLErrors;

p = inputParser;
p.addParameter('behrdir','/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014',@(x) exist(x,'dir'))
p.addParameter('color','w');
p.addParameter('cbrange',[],@(x) length(x) == 2);
p.addParameter('states',1,@isscalar);
p.addParameter('fileprefix','OMI_BEHR_v2-1Arev1_',@isstr);
p.addParameter('clouds','omi',@isstr);
p.addParameter('cloudfraccrit',-1,@isscalar)
p.addParameter('rowanomaly','XTrackFlags',@(x) any(strcmpi(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}))) %Ensure that the rowanomaly value is one of the allowed 4

p.parse(varargin{:});
pout = p.Results;

behr_dir = pout.behrdir;
behr_prefix = pout.fileprefix;
border_col = pout.color;
cbrange = pout.cbrange;
states_bool = pout.states;
cloud_type = pout.clouds;
cloud_frac_crit = pout.cloudfraccrit;
row_anomaly_crit = pout.rowanomaly;

if numel(lonlim) ~= 2 || numel(latlim) ~= 2
    E.badinput('lonlim and latlim must be 2-element vectors');
elseif lonlim(1) >= lonlim(2) || latlim(1) >= latlim(2)
    E.badinput('The lower value must come first in lonlim and latlim')
end

if cloud_frac_crit < 0
    if strcmpi(cloud_type, 'omi')
        cloud_frac_crit = 0.2;
    elseif strcmpi(cloud_type, 'modis')
        cloud_frac_crit = 0;
    else
        E.badinput('Unsupported cloud type: %s', cloud_type)
    end
end
        
    

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

fname = sprintf('%s%04d%02d%02d.mat', behr_prefix, year(date_in), month(date_in), day(date_in));
D = load(fullfile(behr_dir, fname), 'Data');
Data = D.Data;

for d=1:numel(Data)
    Data(d).Areaweight = ones(size(Data(d).BEHRColumnAmountNO2Trop));
    if all(Data(d).Longitude(:) < lonlim(1)) || all(Data(d).Longitude(:) > lonlim(2)) ||...
            all(Data(d).Latitude(:) < latlim(1)) || all(Data(d).Latitude(:) > latlim(2))
        continue
    end
    
    Data(d) = omi_pixel_reject(Data(d), cloud_type, cloud_frac_crit, row_anomaly_crit);
    Data(d).BEHRColumnAmountNO2Trop(Data(d).Areaweight==0) = nan;
    
    figure;
    m_proj('Mercator', 'long', lonlim, 'lat', latlim);
    m_pcolor(Data(d).Longitude, Data(d).Latitude, Data(d).BEHRColumnAmountNO2Trop);
    colorbar;
    if states_bool
        m_states(border_col);
    end
    if ~isempty(cbrange)
        caxis(cbrange);
    end
    m_grid('linestyle','none')
    title(sprintf('Swath %d %s',d,date_in));
end

end

