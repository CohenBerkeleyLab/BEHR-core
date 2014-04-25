%usconic - prepares a map of the US using the Albers-Equal Area Conic
%projection type, drawing the states and coastline in black. Uses the m_map
%package

m_proj('Albers Equal-Area Conic','lon',[-125 -65],'lat',[25 50]);
m_coast('color','k');
m_states('k');
m_grid('linestyle','none');