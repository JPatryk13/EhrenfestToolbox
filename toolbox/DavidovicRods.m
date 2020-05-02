% DAVIDOVICRODS - 
%
%   obj = DavidovicRods(noOfRods, radius, speed), constructor.
%           Input:
%       'noOfRods':
%       'radius':
%       'speed':
%           Output:
%       'obj': object of the class.
%
%   starDisk(obj),
%           Input:
%       'obj':
%
%   gifBoth(obj),
%           Input:
%       'obj':
%
%   Limitations:
%       ...
%
%   Examples:
%       ...
%
%   Use:
%       ...
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, PLOT, GIF, CHANGENOTATIONTYPE
%
%   Patryk Jesionka, 2020
%

classdef DavidovicRods
    properties
        noOfRods {mustBeGreaterThan(noOfRods, 2), mustBeInteger} % Number of rods
        speed = [] % Array of linear speeds of the rotating disc
        radius {mustBeGreaterThan(radius, 0)} % Radius of the inner circle
        
        c = 2.998*10^8; % Speed of light
        
        % Functions describing Lorenz factor and beta angle (described in the constructor)
        gammaFunc = @(v, c) 1./(1-(v.^2./c^2)); % Lorenz factor as a function of speed (v)
        betaFunc = @(g, n) 360./(2*n*g); % Beta angle as a function of Lorenz factor (g)
        
        % Properties described in the constructor
        gamma
        beta
        a
        
        % A function described length of a rod; r - radius, b - beta angle
        rodLenFunc = @(r, b) 2*r*tand(b);
    end
    methods
        function obj = DavidovicRods(noOfRods, radius, speed)
            % User input validation
            obj.noOfRods = noOfRods;
            obj.radius = radius;
            for i = 1:length(speed)
                if ~isnumeric(speed(i))
                    error("Speed of the disk must be a numeric value.")
                end
            end
            if mod(length(speed), 2) == 1
                if ~(length(speed) == 1 || length(speed) == 3)
                    error("Size of the speed array must be 1, 3 or any even number.");
                end
            end
            obj.speed = speed;
            
            obj.gamma = obj.gammaFunc(obj.speed, obj.c); % Lorentz factor(s)
            % Angle(s) between inner (pointing a middle point of a
            % rod) and outer (pointing ends of a rod) circle radii
            obj.beta = obj.betaFunc(obj.gamma, obj.noOfRods);
            
            obj.a = obj.radius./cosd(obj.beta); % Outer circle radius
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
                layout = defineTile(layout, 'title', "v=" + changeNotationType(obj.speed(i), 's') + "m/s, r=" + string(obj.radius) + "m",...
                                            'axesNames', {'x [m]' 'y [m]'},...
                                            'size', axis,...
                                            'legend', 'top-right');
            end
            
            for i = 1:length(obj.speed)
                % Empty array for storing rod's coordinates
                rod = {};
                % Empty array for storing radii pointing the middle of each rod coordinates
                r = {};
                
                % Calculate end-point coordinates for rods and radii
                for j = 1:obj.noOfRods
                    ang = (j-1)*(360/obj.noOfRods);

                    % Calculating coordinates of each of the end points of the rod
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
            
            % Remove legend
            l = findobj('tag','legend');
            delete(l);
            
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
            
            % 
            rodLen = obj.rodLenFunc(obj.radius, obj.beta);
            
            layout = Plot;
            layout = createLayout(layout, 1, 1);
            layout = defineTile(layout, "Rod length vs linear speed", {'Speed [m/s]' 'Rod length [m]'}, [0 0 0 0]);
            layout = addPlot(layout, 1, {speedDomain, rodLenDomain}, '-', 'none', 'r', 0.5, "Change of length of a rod");
            layout = addPlot(layout, 1, {obj.speed, rodLen}, 'none', 'd', 'k', 0.1, "Rod length for an arbitrary" + newline + "chosen speed values");
            drawLayout(layout);
        end
        
        function gifBoth(obj)
            % Code to be finished
        end
    end
end
        