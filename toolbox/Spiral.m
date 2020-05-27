%% SPIRAL
% Generating data required to plot a spiral-shaped trajectory in two- or
% three-dimensional space.
%
%   obj = Spiral(rMin, rMax)
%   obj = Spiral(rMin, rMax, Name, Value) 
%       Constructor, validates user input, generates arrays of 
%       coordinates*, based on the 'fixedCoordinate' decides on the 
%       orientation of the spiral in 3D space and assigns data (with 
%       determined preferred size of the plot) to the structure.
%       *Approach similar to the Circle class'; however, radius range and
%       number of laps is a factor in the coordinates generation.
%
%           Input:
%       'rMin', 'rMax': (required), nonnegative numbers ('rMax' must be
%                       greater than 'rMin'), range of radii in between
%                       which the spiral is plotted.
%       'centrePoint':  [0 0](default), two- or three-dimensional numeric
%                       row array, number of elements determines whether
%                       the plot will be 2D or 3D; it is the centre of a
%                       spiral.
%       'fixedCoordinate': 'z' (default), 'x', 'y' or 'z'; the spiral
%                       plotted in 3D though as a plot has two dimensions -
%                       the letter here states the axis in along which the
%                       circle should be flat. It is related to the
%                       limitation of the class.
%       'noOfLaps':     1 (default), number of laps.
%       'q':            360 (default), must be positive integer, quality 
%                       factor (see Updates).
%           Output:
%       'obj':          object of the class.
%
%   spiral = getSpiral(obj)
%       extract generated data.
%
%           Input:
%       'obj':          object of the class.
%           Output:
%       'spiral':       structure containing 'coordinates' and 'size'
%                       arrays.
%
%   [a, b] = generateSpiral(obj)
%       Private method, performing an actual generation of the spiral
%       coordinates based on the radii range and angle array.
%       
%           Input:
%       'obj':          object of the class
%           Output:
%       [a, b]:         a pair of arrays with coordinates
%
%   Limitations:
%       Generated spiral cannot be oriented freely in 3D space. I.e.
%       the axis going through the centre (perpendicular to the surface the 
%       spiral is plotted on) must be parallel to one of the axes.
%
%   Examples:
%       Functions spiralEx1, spiralEx2 in Examples.m
%
%   Updates:
%       01/03/2020: Added quality factor (constructor). Need for unifying
%           CIRCLE, WAVEFUNCTION and  SPIRAL classes' output number of
%           steps replaced with freedom of choice.
%       23/04/2020: Added input parser.
%
%   Use:
%       Such a structure ('spiral') can be fed into standard MATLAB
%       functions (e.g. plot(), plot3()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, CIRCLE, WAVEFUNCTION, PLOT, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE,
%       FINDLIMITS, CURRENTDENSITY, MAGNETICFLUX
%
%   Patryk Jesionka, 2019
%%

classdef Spiral
    properties
        rMin        % Radius of the inner circle
        rMax        % Radius of the outer circle
        ang         % Angle range of a spiral with defined step
        
       % Structure to be returned with coordinates and default dimensions
       % of the plot
        spiral = struct('coordinates', [], 'size', [])
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
        function obj = Spiral(rMin, rMax, varargin)
            % Define default values
            defaultCentrePoint = [0 0];
            defaultFixedCoordinate = 'z';
            defaultNoOfLaps = 1;
            defaultQ = 360;
            
            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validCentrePoint = @(x) all(isreal(x)) && all(isnumeric(x)) && all(isfinite(x)) && isrow(x) && (eq(length(x), 2) || eq(length(x), 3));
            validFixedCoordinate = @(x) ischar(x) && ismember(x, {'x', 'y', 'z'});
            validNoOfLaps = @(x) ge(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x) && eq(x, floor(x));
            validQ = @(x) validRadius(x) && eq(x, floor(x));
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'rMin', validRadius);
            addRequired(p, 'rMax', validRadius);
            addParameter(p, 'centrePoint', defaultCentrePoint, validCentrePoint);
            addParameter(p, 'fixedCoordinate', defaultFixedCoordinate, validFixedCoordinate);
            addParameter(p, 'noOfLaps', defaultNoOfLaps, validNoOfLaps);
            addParameter(p, 'q', defaultQ, validQ);
            
            parse(p, rMin, rMax, varargin{:});
            
            % Extract variables from the parser
            obj.rMin = p.Results.rMin;
            
            % Check if rMax is greater than rMin
            if p.Results.rMax > obj.rMin
                obj.rMax = p.Results.rMax;
            else
                error("rMax must be greater than rMin.");
            end
            
            centrePoint = p.Results.centrePoint;
            fixedCoordinate = p.Results.fixedCoordinate;
            noOfLaps = p.Results.noOfLaps;
            q = p.Results.q;
            
            % Separate coordinates of spiral's centre
            x = centrePoint(1);
            y = centrePoint(2);
            if length(centrePoint) == 3
                z = centrePoint(3);
            end
            
            % Generating an angle range with step 2*pi/90
            obj.ang = 0:(noOfLaps*2*pi/q):noOfLaps*2*pi;
            
            % Generating coordinate arrays
            [a, b] = obj.generateSpiral;
            c = zeros(1, length(b));
            
            % Decides whether to return 3d or 2d vector to the structure
            if length(centrePoint) == 2
                obj.spiral.coordinates = {x+a y+b};
                
                % Sets optimal size of the plot
                xlim = [x-2*obj.rMax x+2*obj.rMax];
                ylim = [y-2*obj.rMax y+2*obj.rMax];
                obj.spiral.size = [xlim ylim];
            else
                % Based on fixedCoordinate sets orientation of the spiral
                % by assigning a, b and c (c is 'fixed') to appropriate
                % coordinates x, y and z.
                if eq(fixedCoordinate, 'x')
                    obj.spiral.coordinates = {x+c y+a z+b};
                elseif eq(fixedCoordinate, 'y')
                    obj.spiral.coordinates = {x+a y+c z+b};
                else
                    obj.spiral.coordinates = {x+a y+b z+c};
                end
                
                % Sets optimal size of the plot
                xlim = [x-2*obj.rMax x+2*obj.rMax];
                ylim = [y-2*obj.rMax y+2*obj.rMax];
                zlim = [z-2*obj.rMax z+2*obj.rMax];
                obj.spiral.size = [xlim ylim zlim];
            end
        end
        
        function spiral = getSpiral(obj)
            spiral = obj.spiral;
        end
    end
end