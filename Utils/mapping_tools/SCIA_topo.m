%SCIA_topo.m   27 April 2005
%Creates a topography map of the SCIA Mapped location
%Contours are color-coded by GPS elevation using m_tbase.m


    figure, colordef white
    m_proj('mercator', 'long', [-145 -50], 'lat', [20 60]);
    m_coast('patch', [143/255,188/255,143/255], 'edgecolor','none');
    m_grid('box','fancy','tickdir','in'); hold on,
    %m_tbase('contour', [0:500:3000]); colormap(flipud(copper)), colorbar
    %m_gshhs_i('Color','k');
    %plot_states_terrain
    %m_plotbndry('C:\SCIA\mfiles\montana','color','w', 'linewidth', 2);
    %xlabel('Contour Data from UCAR Terrain Base (m)')
    plot_states
    m_gshhs_i('Color','k');