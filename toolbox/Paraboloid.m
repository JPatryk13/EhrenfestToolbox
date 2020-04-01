% PARABOLOID - generating data required to plot a paraboloid/hyperboloid
% surface three-dimensional space.
%
%   obj = Paraboloid(centrePoint, semiAxes, orientation, meshDens), 
%       constructor, validates user input, generates matrices of 
%       coordinates and assigns data (with determined preferred size of the
%       plot) to the structure.
%           Input:
%       'centrePoint': array of numeric values, displacement of the 'tip'
%       of paraboloid relative to the origin of the plot
%       'semiAxes': array of numeric values, length of the semi-minor and
%       semi-major exes of the ellipse. I.e. parameters that dictate the 
%       level of curvature in the x-z and y-z planes respectively.
%       'orientation': '-' or '+' character converted to -1 or 1 integer,
%       determines in which direction curvature opens.
%       'meshDens': positive value, density of the mesh.
%           Output:
%       'obj': object of the class
%
%   paraboloid = getParaboloid(obj), extract generated data
%           Input:
%       'obj': object of the class
%           Output:
%       'paraboloid': structure containing 'coordinates' and 'size' arrays.
%
%   Limitations:
%       Generated paraboloid cannot be oriented freely in 3D space. I.e.
%       The surface can be open upwards or downwards only along the z
%       axis.
%
%   Examples:
%       Plot an elliptical paraboloid opening upwards, starting at the
%       origin.
%           paraboloidHandle = Paraboloid([0 0 0], [1 3], '+', 5);
%           paraboloid = getParaboloid(paraboloidHandle);
%           surf(paraboloid.coordinates{1}, paraboloid.coordinates{2},
%                paraboloid.coordinates{3});
%           axis(paraboloid.size)
%       Plot a hyperboloid with origin at [3 3 1] with a dense grid.
%           paraboloidHandle = Paraboloid([3 3 1], [1 -1], '+', 25);
%           paraboloid = getParaboloid(paraboloidHandle);
%           surf(paraboloid.coordinates{1}, paraboloid.coordinates{2},
%                paraboloid.coordinates{3});
%           axis(paraboloid.size)
%
%   Use:
%       Such a structure ('paraboloid') can be fed into standard MATLAB
%       functions (e.g. surf()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       CIRCLE, SPIRAL, WAVEFUNCTION, PLOT, QUANTUMN, ENERGYAPPROXIMATION,
%       WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE
%
%   Patryk Jesionka, 2019

classdef Paraboloid
    properties
        % x, y, z: coordinates of the paraboloid's centre
        x {mustBeNumeric} = 0
        y {mustBeNumeric} = 0
        z {mustBeNumeric} = 0
        % semi-minor and semi-major axes
        semX {mustBeNumeric} = 1
        semY {mustBeNumeric} = 1
        % Whether it should be directed upwards or downwards
        orientation {mustBeMember(orientation, {'+', '-'})} = '+'
        % How dense the mesh should be
        meshDens {mustBePositive} = 10
        
        % Structure to be returned with coordinates and dimensions of the
        % plot
        paraboloid = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Paraboloid(centrePoint, semiAxes, orientation, meshDens)
            % Input validation
            if ~(length(centrePoint) == 3)
                error("centrePoint array must have 3 entries!");
            end
            obj.x = centrePoint(1);
            obj.y = centrePoint(2);
            obj.z = centrePoint(3);
            if ~(length(semiAxes) == 2)
                error("semiAxes array must have 2 entries");
            end
            obj.semX = semiAxes(1);
            obj.semY = semiAxes(2);
            obj.orientation = orientation;
            obj.meshDens = meshDens;
            
            % Converting orientation to the numerical value (1 or -1)
            numOrnt = str2double(strcat(obj.orientation, '1'));
            
            % Determine semi-major axis
            semiMajor = max([obj.semX obj.semY]);
            
            % "Returns 2-D grid coordinates based on the coordinates contained in
            % vectors x and y. X is a matrix where each row is a copy of x, and Y
            % is a matrix where each column is a copy of y."
            [r, ang] = meshgrid(0:(1/obj.meshDens):semiMajor*10, 0:pi/20:2*pi);

            % Based on the generated r and ang, creates coordinate matrices
            % and decides on the orientation of the surface
            x = r.*cos(ang);
            y = r.*sin(ang);
            z = numOrnt*((x.^2)/obj.semX + (y.^2)/obj.semY);
            
            obj.paraboloid.coordinates = {x+obj.x y+obj.y z+obj.z};
            
            % Generates size of the plot (tile)
            xlim = [obj.x-2*semiMajor obj.x+2*semiMajor];
            ylim = [obj.y-2*semiMajor obj.y+2*semiMajor];
            zlim = sort(numOrnt*[obj.z obj.z+2*semiMajor]);
            
            obj.paraboloid.size = [xlim ylim zlim]; 
        end
        
        function paraboloid = getParaboloid(obj)
            paraboloid = obj.paraboloid;
        end
    end
end

