% CIRCLE - generating data required to plot a circle in two- or
% three-dimensional space.
%
%   obj = Circle(radius)
%   obj = Circle(radius, Name, Value) 
%       constructor, validates user input, generates arrays of coordinates,
%       based on the 'fixedCoordinate' decides on the orientation of the
%       circle in 3D space and assigns data (with determined preferred size 
%       of the plot) to the structure.
%
%           Input parameters' names:
%       'radius':       (required), nonnegative number, radius of a circle 
%                       to create
%       'centrePoint':  [0 0](default), two- or three-dimensional numeric
%                       row array, number of elements determines whether
%                       the plot will be 2D or 3D; it is the centre of a
%                       circle.
%       'fixedCoordinate': 'z' (default), 'x', 'y' or 'z'; the circle
%                       plotted in 3D though as a plot has two dimensions -
%                       the letter here states the axis in along which the
%                       circle should be flat. It is related to the
%                       limitation of the class.
%       'q':            360 (default), must be positive integer, quality 
%                       factor (see Updates).
%           Output:
%       'obj':          object of the class
%
%   circle = getCircle(obj)
%       extract generated data
%
%           Input:
%       'obj':          object of the class
%           Output:
%       'circle':       structure containing 'coordinates' and 'size' 
%                       arrays.
%
%   Limitations:
%       Generated circle cannot be oriented freely in 3D space. I.e.
%       the axis going through the centre (perpendicular to the surface the
%       circle is plotted on) must be parallel to one of the axes.
%
%   Examples:
%       functions circleEx1, circleEx2 in Examples.m
%
%   Updates:
%       01/03/2020: Added quality factor (constructor). Need for unifying
%           CIRCLE, WAVEFUNCTION and  SPIRAL classes' output - number of
%           elements in the coordinates array must be the same when plotted
%           together.
%       22/04/2020: Added input parser.
%
%   Use:
%       Such a structure ('circle') can be fed into standard MATLAB
%       functions (e.g. plot(), plot3()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE
%
%   Patryk Jesionka, 2019
 
classdef Circle
    properties
       % Structure to be returned with coordinates and default dimensions
       % of the plot
       circle = struct('coordinates', [], 'size', [])
    end
    methods
        % Creating an object, validating user input, generating and saving 
        % data in class' properties
        function obj = Circle(radius, varargin)
            % Define default values
            defaultCentrePoint = [0 0];
            defaultFixedCoordinate = 'z';
            defaultQ = 360;
            
            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validCentrePoint = @(x) all(isreal(x)) && all(isnumeric(x)) && all(isfinite(x)) && isrow(x) && (eq(length(x), 2) || eq(length(x), 3));
            validFixedCoordinate = @(x) ischar(x) && ismember(x, {'x', 'y', 'z'});
            validQ = @(x) validRadius(x) && eq(x, floor(x));
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'radius', validRadius);
            addParameter(p, 'centrePoint', defaultCentrePoint, validCentrePoint);
            addParameter(p, 'fixedCoordinate', defaultFixedCoordinate, validFixedCoordinate);
            addParameter(p, 'q', defaultQ, validQ);
            
            parse(p, radius, varargin{:});
            
            % Extract variables from the parser
            radius = p.Results.radius;
            centrePoint = p.Results.centrePoint;
            fixedCoordinate = p.Results.fixedCoordinate;
            q = p.Results.q;
            
            % Separate coordinates of circle's centre
            x = centrePoint(1);
            y = centrePoint(2);
            if length(centrePoint) == 3
                z = centrePoint(3);
            end
 
            % Angle range of a spiral with defined step
            ang = 0:(2*pi/q):2*pi;
            
            % Generate pair of arrays of coordinates
            a = radius*cos(ang);
            b = radius*sin(ang);
            % The third array is fixed
            c = zeros(1, length(b));
 
            % Decides whether to return 3d or 2d vector to the structure
            if length(centrePoint) == 2
                obj.circle.coordinates = {x+a y+b};
                
                % Sets optimal size of the plot
                xlim = findLimits(a, x);
                ylim = findLimits(b, y);
                
                obj.circle.size = [xlim ylim];
            else
                % Based on fixedCoordinate sets orientation of the circle
                % by assigning a, b and c (which is 'fixed') to appropriate
                % coordinates x, y and z.
                if fixedCoordinate == 'x'
                    obj.circle.coordinates = {x+c y+a z+b};
                elseif fixedCoordinate == 'y'
                    obj.circle.coordinates = {x+a y+c z+b};
                else
                    obj.circle.coordinates = {x+a y+b z+c};
                end
                
                % Sets optimal size of the plot
                xlim = findLimits(a, x);
                ylim = findLimits(b, y);
                zlim = findLimits(c, z);
                
                obj.circle.size = [xlim ylim zlim];
            end
        end
        
        % Return structure concerning circle's graph coordinates and its
        % size
        function circle = getCircle(obj)
            circle = obj.circle;
        end
    end
end