% CIRCLE - generating data required to plot a circle in two- or
% three-dimensional space.
%
%   obj = Circle(radius, centrePoint, fixedCoordinate, q), constructor,
%       validates user input, generates arrays of coordinates, based on the
%       'fixedCoordinate' decides on the orientation of the circle in 3D 
%       space and assigns data (with determined preferred size of the plot) 
%       to the structure.
%           Input:
%       'radius': nonnegative number, radius of a circle to create
%       'centrePoint': two- or three-dimensional array of numeric values,
%       dimensionality forces dimensionality of the 'circle' structure's
%       arrays and thus whether the plot will be 2D or 3D; it is the centre 
%       of a circle.
%       'fixedCoordinate': must be a character 'x', 'y' or 'z'; the circle
%       plotted in 3D though as a plot has two dimensions - the letter here
%       states the axis in along which the circle should be flat. It is
%       related to the limitation of the class.
%       'q': must be positive integer, quality factor (see Updates).
%           Output:
%       'obj': object of the class
%
%   circle = getCircle(obj), extract generated data
%           Input:
%       'obj': object of the class
%           Output:
%       'circle': structure containing 'coordinates' and 'size' arrays.
%
%   Limitations:
%       Generated circle cannot be oriented freely in 3D space. I.e.
%       the axis going through the centre (perpendicular to the surface the
%       circle is plotted on) must be parallel to one of the axes.
%
%   Examples:
%       Create a circle in 2D (x-y plane) with radius 5, centred at [0, 0].
%           circleHandle = Circle(5, [0 0], 'z');
%           circle = getCircle(circleHandle);
%           plot(circle.coordinates{1}, circle.coordinates{2});
%           axis(circle.size)
%       Create a circle in 3D (fixed x coordinate) with radius 3.14,
%       centred at [0, 3.14, 3.14] - it is plotted on the y-z surface at
%       x = 0.
%           circleHandle = Circle(3.14, [0 3.14 3.14], 'x');
%           circle = getCircle(circleHandle);
%           plot3(circle.coordinates{1}, circle.coordinates{2}, 
%                circle.coordinates{3});
%           axis(circle.size)
%
%   Updates:
%       01/03/2020: Added quality factor (constructor). Need for unifying
%           CIRCLE, WAVEFUNCTION and  SPIRAL classes' output number of
%           steps replaced with freedom of choice.
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
       % Radius of a circle
       radius {mustBeNonnegative}
       % Constant coordinate - determines orientation of the circle
       fixedCoordinate {mustBeMember(fixedCoordinate, {'x', 'y', 'z'})} = 'z'
       % x, y, z: coordinates of the circle's centre
       x {mustBeNumeric} = 0
       y {mustBeNumeric} = 0
       z {mustBeNumeric} = 0
       
       q {mustBePositive, mustBeInteger} = 360 % Quality factor
 
       % Structure to be returned with coordinates and dimensions of the
       % plot
       circle = struct('coordinates', [], 'size', [])
    end
    methods
        % Creating an object, validating user input, generating and saving 
        % data in class' properties
        function obj = Circle(radius, centrePoint, fixedCoordinate, q)
            % Input validation
            obj.radius = radius;
            obj.fixedCoordinate = fixedCoordinate;
            obj.q = q;
            if ~(length(centrePoint) == 2 || length(centrePoint) == 3)
                error("centrePoint array must have 2 or 3 entries!");
            end
            obj.x = centrePoint(1);
            obj.y = centrePoint(2);
            if length(centrePoint) == 3
                obj.z = centrePoint(3);
            end
 
            % Angle range of a spiral with defined step
            ang = 0:(2*pi/obj.q):2*pi;
            
            % Generate pair of arrays of coordinates
            a = obj.radius*cos(ang);
            b = obj.radius*sin(ang);
            % The third array is fixed
            c = zeros(1, length(b));
 
            % Decides whether to return 3d or 2d vector to the structure
            if length(centrePoint) == 2
                obj.circle.coordinates = {obj.x+a obj.y+b};
                % Sets optimal size of the plot
                xlim = [obj.x-2*obj.radius obj.x+2*obj.radius];
                ylim = [obj.y-2*obj.radius obj.y+2*obj.radius];
                obj.circle.size = [xlim ylim];
            else
                % Based on fixedCoordinate sets orientation of the circle
                % by assigning a, b and c (which is 'fixed') to appropriate
                % coordinates x, y and z.
                if obj.fixedCoordinate == 'x'
                    obj.circle.coordinates = {obj.x+c obj.y+a obj.z+b};
                elseif obj.fixedCoordinate == 'y'
                    obj.circle.coordinates = {obj.x+a obj.y+c obj.z+b};
                else
                    obj.circle.coordinates = {obj.x+a obj.y+b obj.z+c};
                end
                % Sets optimal size of the plot
                xlim = [obj.x-2*obj.radius obj.x+2*obj.radius];
                ylim = [obj.y-2*obj.radius obj.y+2*obj.radius];
                zlim = [obj.z-2*obj.radius obj.z+2*obj.radius];
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