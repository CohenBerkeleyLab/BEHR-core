%function [ pix_lon, pix_lat, pix_loncorn, pix_latcorn ] = tempo_generator(  )
function [ Data ] = tempo_generator(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTANTS %%%%%
%%%%%%%%%%%%%%%%%%%%%

r_E = 6375; % radius earth in km
tempo_alt = 35786; % tempo altitude above earth's surface in km
tempo_lon = -100;
tempo_lat = 0;
for_centlon = -100; % center coordinates for the field of regard
for_centlat = 36.5;
ns_angle = 0.0036 / 180 * pi; % N/S field of view angle for one pixel
ew_angle = 0.0070 / 180 * pi; % E/W arc subtended for one pixel during scan.

az_0 = 80 / 180 * pi;
elev_0 = 5.82 / 180 * pi;

tempo_xyz = zeros(3,1);
[tempo_xyz(1), tempo_xyz(2), tempo_xyz(3)] = sph2cart(deg2rad(tempo_lon), deg2rad(tempo_lat), tempo_alt + r_E);

n_ew = 3000;
n_ns = 2500;



grid_loncorn = 1000*ones(n_ns, n_ew);
grid_latcorn = 1000*ones(n_ns, n_ew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRESSIVE VARS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for a=1:(n_ns/2)
    if mod(a,100) == 0; fprintf('a = %d of %d. Time elapsed = %f\n',a,n_ns/2,toc); end
    for b=1:(n_ew/2)
        %if mod(b,100)==0; fprintf('\t b = %d. Time = %f\n',b,toc); end
        i = n_ns/2 + a;
        j = n_ew/2 + b;
        [grid_loncorn(i,j), grid_latcorn(i,j)] = calc_pix_corners(a,b);
        
        i = n_ns/2 + 1 - a;
        j = n_ew/2 + 1 - b;
        [grid_loncorn(i,j), grid_latcorn(i,j)] = calc_pix_corners(-a,-b);
        
        i = n_ns/2 + 1 - a;
        j = n_ew/2 + b;
        [grid_loncorn(i,j), grid_latcorn(i,j)] = calc_pix_corners(-a,b);
        
        i = n_ns/2 + a;
        j = n_ew/2 + 1 - b;
        [grid_loncorn(i,j), grid_latcorn(i,j)] = calc_pix_corners(a,-b);
    end
end

% Now convert the grid into pixel lat/lon centers and corners.
sz = size(grid_loncorn)-1;
Data.Longitude = nan(sz);
Data.Latitude = nan(sz);
Data.Loncorn = nan(4, sz(1), sz(2));
Data.Latcorn = nan(4, sz(1), sz(2));
for a=1:sz(1)
    if mod(a,100) == 0; fprintf('a = %d of %d for pixel writing. Time elapsed = %f\n',a,sz(1),toc); end
    for b=1:sz(2)
        lons = reshape(grid_loncorn(a:a+1,b:b+1),4,1);
        lats = reshape(grid_latcorn(a:a+1,b:b+1),4,1);
        Data.Longitude(a,b) = mean(lons);
        Data.Latitude(a,b) = mean(lats);
        if ~isnan(mean(lons))
            Data.Loncorn(:,a,b) = lons;
            Data.Latcorn(:,a,b) = lats;
        end
    end
end


    function [loncorn, latcorn] = calc_pix_corners(a,b)
        ew_angle_i = (b - 0.5*sign(b)) * ew_angle + az_0;
        ns_angle_j = (a - 0.5*sign(a)) * ns_angle + elev_0;
        %ew_angle_i = az_0;
        %ns_angle_j = elev_0;
        
        % The position of the grid points is found by solving for the
        % intersection of a ray propagating from a pixel corner of TEMPO to
        % the Earth surface. The intersection is defined as the point where
        % the (x,y,z) coordinates from that ray are a solution to the
        % equation of a sphere with the radius of the earth. The ray is
        % defined by three parametric equations in r, where r is the
        % distance from TEMPO.
        %
        % System of equations:
        %   x^2 + y^2 + z^2 = r_E^2
        %   x = x_t + m_x * r
        %   y = y_t + m_y * r
        %   z = z_t + m_z * r
        %
        % x_t, y_t, and z_t are the cartesian coordinates for TEMPO in a
        % coordinate system with the center of the earth at the origin. The
        % rays are defined above relative to the line of sight from TEMPO
        % to the center of its domain.
        
        [m.x, m.y, m.z] = sph2cart(ew_angle_i, ns_angle_j, 1);
        a_q = m.x^2 + m.y^2 + m.z^2;
        b_q = 2*(m.x*tempo_xyz(1) + m.y*tempo_xyz(2) + m.z*tempo_xyz(3));
        c_q = sum(tempo_xyz.^2) - r_E^2;
        r = (-b_q - sqrt( b_q^2 - 4*a_q*c_q) ) / 2*a_q;
        if r < 0
            error('quad_form:solv','r should not be < 0')
        elseif imag(r)
            loncorn = nan;
            latcorn = nan;
            return
        end
        
        x = tempo_xyz(1) + m.x * r;
        y = tempo_xyz(2) + m.y * r;
        z = tempo_xyz(3) + m.z * r;
        
        [loncorn, latcorn] = cart2sph(x, y, z);
        loncorn = loncorn / pi * 180;
        latcorn = latcorn / pi * 180;
        
%         r = 0.99*tempo_alt;
%         r_xyz = zeros(3,1);
%         last_norm = inf;
%         while true
%             r = r+1;
%             [r_xyz(1), r_xyz(2), r_xyz(3)] = sph2cart(ew_angle_i, ns_angle_j, r);
%             r_xyz = r_xyz + tempo_xyz;
%             if norm(r_xyz) < r_E
%                 break
%             elseif norm(r_xyz) > last_norm
%                 error('huh:weird','r is increasing')
%             end
%             last_norm = norm(r_xyz);
%             last_rxyz = r_xyz;
%             last_r = r;
%         end
%         while true
%             r = r-0.01;
%             [r_xyz(1), r_xyz(2), r_xyz(3)] = sph2cart(ew_angle_i, ns_angle_j, r);
%             r_xyz = r_xyz + tempo_xyz;
%             if norm(r_xyz) > r_E
%                 break
%             end
%         end
%         [loncorn, latcorn] = cart2sph(r_xyz(1), r_xyz(2), r_xyz(3));
%         loncorn = loncorn / pi * 180;
%         latcorn = latcorn / pi * 180;
    end

end

