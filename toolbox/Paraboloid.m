%% PARABOLOID
% generating data required to plot a paraboloid/hyperboloid surface
% three-dimensional space.
%
%   obj = Paraboloid(semiAxes)
%   obj = Paraboloid(semiAxes, Name, Value)
%       constructor, validates user input, generates matrices of 
%       coordinates and assigns data (with determined size of the plot)
%       to the structure.
%
%           Input:
%       'semiAxes':     (required), array of numeric values, length of the
%                       semi-minor and semi-major exes of the ellipse,
%                       i.e. parameters that dictate the level of curvature
%                       in the x-z and y-z planes respectively.
%       'centrePoint':  [0 0 0] (default), array of numeric values, 
%                       displacement of the 'tip' of paraboloid relative
%                       to the origin of the plot.
%       'orientation':  '+' (default), '-' or '+' character converted to
%                       -1 or 1 integer, determines in which direction
%                       curvature opens.
%       'type':         'normal' (default), 'normal' or 'hyperbolic' allows
%                       for plotting a hyperbolic paraboloid as well.
%       'meshDens':     1 (default), positive value, density of the mesh.
%       'rLim':         0 (default), if specified different than zero and
%                       smaller than semiMajor then plots the paraboloid
%                       only to the distance from the peak specified by
%                       rLim.
%           Output:
%       'obj':          object of the class
%
%   paraboloid = getParaboloid(obj)
%       extract generated data.
%
%           Input:
%       'obj':          object of the class
%           Output:
%       'paraboloid':   structure containing 'coordinates' and 'size'
%                       arrays.
%
%   Limitations:
%       Generated paraboloid cannot be oriented freely in 3D space, i.e.
%       The surface can be open upwards or downwards only along the z
%       axis.
%
%   Examples:
%       Functions paraboloidEx1, paraboloidEx2 from Examples.m
%
%   Updates:
%       23/04/2020: Added input parser.
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
%%

classdef Paraboloid
    properties
        % Structure to be returned with coordinates and dimensions of the
        % plot
        paraboloid = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Paraboloid(semiAxes, varargin)
            % Define default values
            defaultCentrePoint = [0 0 0];
            defaultOrientation = '+';
            defaultType = 'normal';
            defaultMeshDens = 10;
            defaultRLim = 0;
            
            % Validation functions
            validSemiAxes = @(x) all(gt(x, 0)) && all(isreal(x)) && all(isnumeric(x)) && all(isfinite(x)) && isrow(x) && eq(length(x), 2);
            validCentrePoint = @(x) all(isreal(x)) && all(isnumeric(x)) && all(isfinite(x)) && isrow(x) && eq(length(x), 3);
            validOrientation = @(x) ismember(x, {'+', '-'});
            validType = @(x) ismember(x, {'normal', 'hyperbolic'});
            validMeshDens = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validRLim = @(x) ge(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'semiAxes', validSemiAxes);
            addParameter(p, 'centrePoint', defaultCentrePoint, validCentrePoint);
            addParameter(p, 'orientation', defaultOrientation, validOrientation);
            addParameter(p, 'type', defaultType, validType);
            addParameter(p, 'meshDens', defaultMeshDens, validMeshDens);
            addParameter(p, 'rLim', defaultRLim, validRLim);
            
            parse(p, semiAxes, varargin{:});
            
            % Extract variables from the parser
            semiAxes = p.Results.semiAxes;
            centrePoint = p.Results.centrePoint;
            orientation = p.Results.orientation;
            type = p.Results.type;
            meshDens = p.Results.meshDens;
            rLim = p.Results.rLim;
            
            % Separate coordinates of paraboloid's centre and lengths of
            % semi axes
            x = centrePoint(1);
            y = centrePoint(2);
            z = centrePoint(3);
            semX = semiAxes(1);
            semY = semiAxes(2);
            
            % Converting orientation to the numerical value (1 or -1)
            numOrnt = str2double(strcat(orientation, '1'));
            
            % Change sign of one of the semi axes when appropriate
            if length(type) == 10
                semY = 0 - semY;
            end
            
            % Determine semi-major axis
            semiMajor = max([semX semY]);
            
            % "Returns 2-D grid coordinates based on the coordinates contained in
            % vectors x and y. X is a matrix where each row is a copy of x, and Y
            % is a matrix where each column is a copy of y."
            [r, ang] = meshgrid(0:(1/meshDens):semiMajor, 0:pi/20:2*pi);
            
            % Limit the displayed area of the paraboloid if rLim is
            % different than 0 and smaller than the semi major axis
            if ~eq(rLim, 0) && gt(semiMajor, rLim)
                for i = length(r(1,:)):-1:1
                    % Loop from the last to the first element of the radius
                    % array removing columns with values greater than
                    % specified radius.
                    last = length(r(1,:));
                    if r(1,last) <= rLim
                        % Break the loop if the radius value is equal or
                        % lower then the radius specified
                        break
                    else
                        % Remove the last column from the radius array and
                        % the angle array to keep the same size
                        r(:,last) = [];
                        ang(:,last) = [];
                    end
                end
            end
            
            % Based on the generated r and ang, creates coordinate matrices
            % and decides on the orientation of the surface
            a = r.*cos(ang);
            b = r.*sin(ang);
            c = numOrnt*((a.^2)/semX + (b.^2)/semY);
            
            obj.paraboloid.coordinates = {a+x b+y c+z};
            
            % Generates size of the plot (tile)
            xlim = findLimits(a, x);
            ylim = findLimits(b, y);
            zlim = findLimits(c, z);
            
            obj.paraboloid.size = [xlim ylim zlim]; 
        end
        
        function paraboloid = getParaboloid(obj)
            paraboloid = obj.paraboloid;
        end
    end
end

