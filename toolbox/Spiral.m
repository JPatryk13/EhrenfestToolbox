% SPIRAL - generating data required to plot a spiral-shaped trajectory
% in two- or three-dimensional space.
%
%   obj = Spiral(rMin, rMax, centrePoint, fixedCoordinate, noOfLaps), 
%       constructor, validates user input, generates arrays of 
%       coordinates*, based on the 'fixedCoordinate' decides on the 
%       orientation of the spiral in 3D space and assigns data (with 
%       determined preferred size of the plot) to the structure.
%       *Approach similar to the Circle class'; however, radius range and
%       number of laps is a factor in the coordinates generation.
%           Input:
%       'rMin', 'rMax': nonnegative numbers ('rMax' must be greater than
%       'rMin'), range of radii in between which the spiral is plotted.
%       'centrePoint': two- or three-dimensional array of numeric values,
%       dimensionality forces dimensionality of the 'spiral' structure's
%       arrays and thus whether the plot will be 2D or 3D; it is the entre 
%       of a spiral.
%       'fixedCoordinate': must be a character 'x', 'y' or 'z'; the spiral
%       plotted in 3D though as a plot has two dimensions - the letter here
%       states the axis in along which the spiral should be flat. It is
%       related to the limitation of the class.
%           Output:
%       'obj': object of the class
%
%   spiral = getSpiral(obj), extract generated data
%           Input:
%       'obj': object of the class
%           Output:
%       'spiral': structure containing 'coordinates' and 'size' arrays.
%
%   [a, b] = generateSpiral(obj), private method, performing an actual
%       generation of the spiral coordinates based on the radii range and 
%       angle array.
%           Input:
%       'obj': object of the class
%           Output:
%       [a, b]: a pair of arrays with coordinates
%
%   Limitations:
%       Generated spiral cannot be oriented freely in 3D space. I.e.
%       the axis going through the centre (perpendicular to the surface the 
%       spiral is plotted on) must be parallel to one of the axes.
%
%   Examples:
%       Create a spiral in 2D (x-y plane) with radii range from 3 to 5, 
%       centred at [0, 0] and 1 lap.
%           spiralHandle = Spiral(3, 5, [0 0], 'z', 1);
%           spiral = getSpiral(spiralHandle);
%           plot(spiral.coordinates{1}, spiral.coordinates{2});
%           axis(spiral.size)
%       Create a spiral in 3D (fixed x coordinate) with radius ranging from
%       6.28 to 12.56 centred at [0, 3.14, 3.14] - it is plotted on the y-z 
%       surface at x = 0 - and 4 laps.
%           spiralHandles = Spiral(6.28, 12.56, [0 3.14 3.14], 'x', 4);
%           spiral = getSpiral(spiralHandle);
%           plot3(spiral.coordinates{1}, spiral.coordinates{2}, 
%                 spiral.coordinates{3});
%           axis(spiral.size)
%
%   Use:
%       Such a structure ('spiral') can be fed into standard MATLAB
%       functions (e.g. plot(), plot3()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, CIRCLE, WAVEFUNCTION, PLOT, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE
%
%   Patryk Jesionka, 2019

classdef Spiral
    properties
        % rMax (later r), rMin: 'starting' and 'ending' radii of the spiral
        rMin {mustBePositive}
        rMax {mustBePositive}
        % Constant coordinate - determines orientation of the spiral
        fixedCoordinate {mustBeMember(fixedCoordinate, {'x', 'y', 'z'})} = 'z'
        % laps: number of full rotations
        noOfLaps {mustBeNonnegative}
        % x, y, z: coordinates of the spiral's centre
        x {mustBeNumeric} = 0
        y {mustBeNumeric} = 0
        z {mustBeNumeric} = 0
        % Structure to be returned with coordinates and dimensions of the
        % plot
        spiral = struct('coordinates', [], 'size', [])
        % angle range of a spiral with defined step
        ang
    end
    methods (Access = private)
        function [a, b] = generateSpiral(obj)
            % Function returns a pair of arrays with coordinates

            % Array of radii going from rMax to rMin with a step dr
            % Step dr normalised so that the r and ang arrays both have the
            % same size
            dr = (obj.rMax - obj.rMin)/(length(obj.ang)-1);
            r = obj.rMax:-dr:obj.rMin;
            
            % Populating arrays with coordinates
            a = r.*cos(obj.ang);
            b = r.*sin(obj.ang);
        end
    end
    methods (Access = public)
        function obj = Spiral(rMin, rMax, centrePoint, fixedCoordinate, noOfLaps)
            % Input validation
            obj.rMin = rMin;
            obj.rMax = rMax;
            if ~(obj.rMin < obj.rMax)
                error("rMin must be smaller than rMax!");
            end
            obj.fixedCoordinate = fixedCoordinate;
            obj.noOfLaps = noOfLaps;
            if ~(length(centrePoint) == 2 || length(centrePoint) == 3)
                error("centrePoint array must have 2 or 3 entries!");
            end
            obj.x = centrePoint(1);
            obj.y = centrePoint(2);
            if length(centrePoint) == 3
                obj.z = centrePoint(3);
            end
            
            % Generating an angle range with step 2*pi/90
            obj.ang = 0:(2*pi/90):obj.noOfLaps*2*pi;
            
            % Generating coordinate arrays
            [a, b] = obj.generateSpiral;
            c = zeros(1, length(b));
            
            % Decides whether to return 3d or 2d vector to the structure
            if length(centrePoint) == 2
                obj.spiral.coordinates = {obj.x+a obj.y+b};
                % Sets optimal size of the plot
                xlim = [obj.x-2*obj.rMax obj.x+2*obj.rMax];
                ylim = [obj.y-2*obj.rMax obj.y+2*obj.rMax];
                obj.spiral.size = [xlim ylim];
            else
                % Based on fixedCoordinate sets orientation of the spiral
                % by assigning a, b and c (c is 'fixed') to appropriate
                % coordinates x, y and z.
                if obj.fixedCoordinate == 'x'
                    obj.spiral.coordinates = {obj.x+c obj.y+a obj.z+b};
                elseif obj.fixedCoordinate == 'y'
                    obj.spiral.coordinates = {obj.x+a obj.y+c obj.z+b};
                else
                    obj.spiral.coordinates = {obj.x+a obj.y+b obj.z+c};
                end
                % Sets optimal size of the plot
                xlim = [obj.x-2*obj.rMax obj.x+2*obj.rMax];
                ylim = [obj.y-2*obj.rMax obj.y+2*obj.rMax];
                zlim = [obj.z-2*obj.rMax obj.z+2*obj.rMax];
                obj.spiral.size = [xlim ylim zlim];
            end
        end
        
        function spiral = getSpiral(obj)
            spiral = obj.spiral;
        end
    end
end