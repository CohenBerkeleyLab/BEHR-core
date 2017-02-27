classdef NO2WindSectors < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        MinWindToAll = true;
        MinWindCrit = 0.1;
        Grid
    end
    properties(SetAccess = protected)
        E
        NE
        N
        NW
        W
        SW
        S
        SE
    end
    properties(Hidden, SetAccess = protected)
        iE = 1;
        iNE = 1;
        iN = 1;
        iNW = 1;
        iW = 1;
        iSW = 1;
        iS = 1;
        iSE = 1;
    end
    properties(Constant)
        directions = {'W','SW','S','SE','E','NE','N','NW'};
        theta_bin_centers = [-180, -135, -90, -45, 0, 45, 90, 135];
    end
    properties(Hidden, Access = protected, Constant)
        theta_bin_edges = [-180, -157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5];
    end
    
    methods
        function obj = NO2WindSectors(Grid, ntimes, minwindtoall)
            if ~isa(Grid,'GlobeGrid')
                error('NO2WindSectors:bad_input','GRID must be an instance of GlobeGrid')
            elseif ~isscalar(ntimes) || ~isnumeric(ntimes) || ntimes < 1 || mod(ntimes,1) ~= 0
                error('NO2WindSectors:bad_input','NTIMES must be a scalar whole number >= 1')
            elseif exist('minwindtoall','var') && (~isscalar(minwindtoall) || ~islogical(minwindtoall))
                error('NO2WindSectors:bad_input','MINWINDTOALL must be a scalar logical (if given)')
            end
            obj.Grid = Grid;
            sz = size(Grid.GridLon);
            init_grid = nan([sz, ntimes]);
            for a=1:numel(NO2WindSectors.directions)
                obj.(NO2WindSectors.directions{a}) = init_grid;
            end
            if exist('minwindtoall','var')
                obj.MinWindToAll = minwindtoall;
            end
        end
        
        %%%%%%%%%%%
        % Setters %
        %%%%%%%%%%%
        function obj = set.MinWindToAll(obj, val) %#ok<MCHV2>
            if ~isscalar(val) || ~islogical(val)
                error('NO2WindSectors:bad_value','MinWindToAll must be a scalar boolean value')
            end
            obj.MinWindToAll = val;
        end
        
        function obj = set.MinWindCrit(obj, val)
            if ~isscalar(val) || ~isnumeric(val) || val < 0
                error('NO2WindSectors:bad_value', 'MinWindCrit must be a scalar number >= 0')
            end
            obj.MinWindCrit = val;
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Primary methods %
        %%%%%%%%%%%%%%%%%%%
        function AddDataToDirection(obj, data, theta, windvel)
            if obj.MinWindToAll && ~exist('windvel','var')
                error('NO2WindSectors:bad_input','If the MinWindToAll boolean is true, must provide a wind speed')
            elseif ~isscalar(theta) || ~isnumeric(theta) || theta < -180 || theta > 180
                error('NO2WindSectors:bad_input','THETA must be a scalar number between -180 and +180')
            elseif obj.MinWindToAll 
                if ~isscalar(windvel) || ~isnumeric(windvel)
                    error('NO2WindSectors:bad_input','WINDVEL must be a scalar number')
                elseif windvel < 0
                    warning('Negative wind speed given')
                end
            end
            
            if obj.MinWindToAll && windvel < obj.MinWindCrit
                for a=1:numel(NO2WindSectors.directions)
                    obj.add_to_direction(data, NO2WindSectors.directions{a});
                end
            else
                % Find which bin this belongs in
                if theta < NO2WindSectors.theta_bin_edges(2) || theta > NO2WindSectors.theta_bin_edges(end)
                    bin = 'W';
                else
                    theta_xx = NO2WindSectors.theta_bin_edges(1:end-1) <= theta & NO2WindSectors.theta_bin_edges(2:end) > theta;
                    bin = NO2WindSectors.directions{theta_xx};
                end
                obj.add_to_direction(data, bin);
            end
        end
        
        function TrimSectorArrays(obj)
            % Removes time slices of the arrays that are all NaNs, but
            % leaves at least one regardless
            for a=1:numel(NO2WindSectors.directions)
                dbin = NO2WindSectors.directions{a};
                sz = size(obj.(dbin));
                nans = all(isnan(reshape(obj.(dbin), sz(1)*sz(2),[])),1);
                if all(nans)
                    nans(1) = false; % leave at least one
                end
                obj.(dbin)(:,:,nans) = [];
            end
        end
    end
    
    methods(Hidden, Access = protected)
        function add_to_direction(obj, data, direction)
            sz = size(obj.(direction));
            if ~isequal(size(data), sz(1:2))
                error('NO2WindSectors:bad_input','DATA is a different size than the sector arrays were initialized to')
            end
            dirindfield = ['i', direction];
            i = obj.(dirindfield);
            obj.(direction)(:,:,i) = data;
            obj.(dirindfield) = i+1;
        end
    end
    
end

