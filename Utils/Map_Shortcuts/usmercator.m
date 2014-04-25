%usmercator - prepares a map of the US using the Mercator
%projection type, drawing the states and coastline in black. Uses the m_map
%package.  Lat/lon boundaries extended; useful for visualizing fuller
%satellite swaths

m_proj('Mercator','lonj',[-125 -65],'lat',[25 50]);
m_coast('color','k');
m_states('k');
m_grid('linestyle','none');