function [ b ] = PlotSwathBoundaries( Data )
%Plots the boundaries of the OMI swaths in "Data" on a map of the US
%   Josh Laughner 20 Mar 2014 <joshlaugh5@gmail.com>

%lines = ['-ro','--bs',':gv','-.md','-ch'];
lines = ['-', '--', ':', '-.'];
linecolors = ['r', 'b', 'g', 'm'];
markers = ['o', 's', 'v', 'd'];

b=length(Data);
 lats = zeros(b,4); lons = zeros(b,4);
for d=1:b
    r1 = find(Data(d).Row == min(Data(d).Row));
    r2 = find(Data(d).Row == max(Data(d).Row));
    
    corners = [min(r1), max(r1), min(r2), max(r2)];
    for c = 1:4
        lats(b,c) = Data(d).Latitude(corners(c));
        lons(b,c) = Data(d).Longitude(corners(c));
    end
end

m_proj('Albers Equal-Area Conic','lon',[-125 -65],'lat',[25 50]);
m_coast;
m_states('b');

for a=1:b
    fprintf('Plotting %u\n',a);
    m_box(lons(a,:),lats(a,:),'mode','corners','linecolor',linecolors(mod(a,4)+1),'marker',markers(mod(a,4)+1),'linestyle',lines(mod(a,4)+1),'linewidth',3);
    pause
end

end

