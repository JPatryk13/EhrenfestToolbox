%% DAVIDOVICRODS
% Plotting meter rods tangent to the circle with defined radius. Input
% parameter, speed defines how much each rod should shrink. The result -
% plot/gif - shows the Lorenz contraction phenomenon affecting rigid body
% rotating with relativistic speed. 
%
%   obj = DavidovicRods(radius, speed)
%   obj = DavidovicRods(radius, speed, 'noOfRods', noOfRodsVal)
%       Constructor, validates input and performs initial calculations
%       carried out throughout the class' functions.
%
%           Input:
%       'radius':       (required), radius of the star disc.
%       'speed':        (required), linear speed of the rotating disc.
%                       Can accept an array of values if many plots are
%                       required to show the change across the speed range.
%       'noOfRods':     10 (default), allows to customise plotted star
%                       disc.
%           Output:
%       'obj':          object of the class.
%
%   starDisk(obj)
%       Function plots graph(s) of star disc with respective rods' lengths.
%       Uses PLOT class to designate tiled layout.
%
%           Input:
%       'obj':          object of the class.
%
%   rodLength(obj)
%       Funtion plots linear regression of the rod length across the whole
%       allowed speed range. Adds markers for user-defined speed values.
%
%           Input:
%       'obj':          object of the class.
%
%   gifBoth(obj)
%       Function creates an animation of the changing rod length w.r.t
%       speed against the rods' length regression - two-plot tiled layout.
%       >> Function to be created <<
%
%           Input:
%       'obj':          object of the class.
%
%   Limitations:
%       N/A
%
%   Examples:
%       Function davidovicRodsEx from Examples.m
%
%   Updates:
%       03/05/2020: Added input parser.
%
%   Use:
%       Displaying Lorenz contraction effect on the circumference of a
%       rigid body.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, PLOT, GIF, CHANGENOTATIONTYPE,
%       FINDLIMITS, CURRENTDENSITY, MAGNETICFLUX
%
%   Patryk Jesionka, 2020
%%

classdef DavidovicRods
    properties
        noOfRods            % Number of rods
        
        speed = []          % Array of linear speeds of the rotating disc
        radius              % Radius of the inner circle
        
        c = 2.998*10^8;     % Speed of light
        
        % Functions describing Lorenz factor and beta angle (described in
        % the constructor).
        gammaFunc = @(v, c) 1./(1-(v.^2./c^2));
        betaFunc = @(g, n) 360./(2*n*g);
        
        % Properties described in the constructor
        gamma
        beta
        a
        
        % A function described length of a rod; r - radius, b - beta angle
        rodLenFunc = @(r, b) 2*r*tand(b);
    end
    methods
        function obj = DavidovicRods(radius, speed, varargin)
            % Define default values
            defaultNoOfRods = 10;
            
            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validSpeed = @(x) all(ge(x, 0)) && all(ge(3*10^8, x)) &&  all(isreal(x)) && all(isnumeric(x)) &&...
                              (isscalar(x) || (isrow(x) && ~iscell(x) && (ismember(length(x), [1 3]) || eq(mod(length(x), 2), 0))));
            validNoOfRods = @(x) gt(x, 2) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x) && eq(floor(x), x);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'radius', validRadius);
            addRequired(p, 'speed', validSpeed);
            addOptional(p, 'noOfRods', defaultNoOfRods, validNoOfRods);
            
            parse(p, radius, speed, varargin{:});
            
            % Extract variables from the parser
            obj.radius = p.Results.radius;
            obj.speed = p.Results.speed;
            obj.noOfRods = p.Results.noOfRods;
            
            % Lorentz factor(s)
            obj.gamma = obj.gammaFunc(obj.speed, obj.c);
            
            % Angle(s) between inner (pointing a middle point of a
            % rod) and outer (pointing ends of a rod) circle radii
            obj.beta = obj.betaFunc(obj.gamma, obj.noOfRods);
            
            % Outer circle radius
            obj.a = obj.radius./cosd(obj.beta);
        end
        
        function starDisk(obj)
            % Obtain size of the layout
            if length(obj.speed) == 1
                m = 1;
                n = 1;
            elseif length(obj.speed) == 3
                m = 1;
                n = 3;
            else
                m = length(obj.speed)/2;
                n = 2;
            end
            
            % Length of axes
            axis = [-2*obj.radius 2*obj.radius -2*obj.radius 2*obj.radius];
            
            % Create layout
            layout = Plot;
            layout = createLayout(layout, m, n);
            for i = 1:length(obj.speed)
                title = "v=" + changeNotationType(obj.speed(i), 's') + "m/s, r=" + string(obj.radius) + "m";
                layout = defineTile(layout, 'title', title,...
                                            'axesNames', {'x [m]' 'y [m]'},...
                                            'size', axis,...
                                            'legend', 'none');
            end
            
            for i = 1:length(obj.speed)
                % Empty array for storing rod's coordinates
                rod = {};
                % Empty array for storing radii pointing the middle of each
                % rod coordinates
                r = {};
                
                % Calculate end-point coordinates for rods and radii
                for j = 1:obj.noOfRods
                    ang = (j-1)*(360/obj.noOfRods);

                    % Calculating coordinates of each of the end points of
                    % the rod
                    x1 = obj.a(i)*sind(ang - obj.beta(i));
                    y1 = obj.a(i)*cosd(ang - obj.beta(i));
                    x2 = obj.a(i)*sind(ang + obj.beta(i));
                    y2 = obj.a(i)*cosd(ang + obj.beta(i));

                    % Calculating coordinates of middle points of the rod
                    x0 = (x1 + x2)/2;
                    y0 = (y1 + y2)/2;

                    % Generating coordinate array for a rod
                    rod(j,:) = {[x1 x2] [y1 y2]};
                    r(j,:) = {[0 x0] [0 y0]};
                end
                
                % Add plots to a tile
                for j = 1:obj.noOfRods
                    layout = addPlot(layout, i, rod(j,:), 'lineSpec', '-k',...
                                                          'lineWidth', 4,...
                                                          'name', "Tangent");
                    layout = addPlot(layout, i, r(j,:), 'lineSpec', '-.b',...
                                                        'name', "Radius");
                end
            end
            
            % Draw layout
            drawLayout(layout);
            
            % Add title to the layout
            t = findobj('type', 'tiledlayout');
            title(t, "Davidovic rods");
        end
        
        function rodLength(obj)
            % Difference between each speed value required to obtain
            % 360-element array of speed values ranging between 0 and
            % 0.99c
            deltaSpeed = obj.c/359;
            
            % Domains (range of values) of speed, Lorenz factor and beta
            % angle - nested function separated for the sake of readability
            speedDomain = 0:deltaSpeed:obj.c;
            gammaDomain = obj.gammaFunc(speedDomain, obj.c);
            betaDomain = obj.betaFunc(gammaDomain, obj.noOfRods);
            
            % Rod length changing from maximum down to 0
            rodLenDomain = obj.rodLenFunc(obj.radius, betaDomain);
            
            rodLen = obj.rodLenFunc(obj.radius, obj.beta);
            
            layout = Plot;
            layout = createLayout(layout, 1, 1);
            layout = defineTile(layout, 'title', "Rod length vs linear speed",...
                                        'axesNames', {'Speed [m/s]'...
                                        'Rod length [m]'});
            % Linear regression of the rod length
            layout = addPlot(layout, 1, {speedDomain, rodLenDomain}, 'color', 'r',...
                                                                     'name', "Change of length of a rod");
            % Markers for each defined speed step
            layout = addPlot(layout, 1, {obj.speed, rodLen}, 'lineSpec', 'o',...
                                                             'markerColor', 'k',...
                                                             'name', "Rod length for" + newline + "chosen speed values");
            % Vertical and horizontal lines pointing markers
            for i = 1:length(obj.speed)
                layout = addPlot(layout, 1, {[0 obj.speed(i)] [rodLen(i) rodLen(i)]}, 'lineSpec', '--k');
                layout = addPlot(layout, 1, {[obj.speed(i) obj.speed(i)] [0 rodLen(i)]}, 'lineSpec', '--k');
            end
            drawLayout(layout);
        end
        
        function gifBoth(obj)
            % Code to be finished
        end
    end
end
        