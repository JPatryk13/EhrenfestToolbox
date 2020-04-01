% DAVIDOVICRODS - 
%
%   obj = Gif, constructor.
%           Output:
%       'obj': object of the class.
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
%       ENERGYAPPROXIMATION, WAVE, PLOT, GIF
%
%   Patryk Jesionka, 2020
%

classdef DavidovicRods
    properties
        noOfRods {mustBeGreaterThan(noOfRods, 2), mustBeInteger} % Number of rods
        speed = [] % Array of linear speeds of the rotating disc
        radius {mustBeGreaterThan(radius, 0)} % Radius of the inner circle
        c = 2.998*10^8; % Speed of light
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
            
            % Circle
            circleHandle = Circle(radius, [0 0], 'z', 360);
            circle = getCircle(circleHandle);
            
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
            
            % Create layout
            layout = Plot;
            layout = createLayout(layout, m, n);
            for i = 1:length(obj.speed)
                layout = defineTile(layout, "v=" + changeNotationType(obj.speed(i), 's') + "m/s, r=" + string(obj.radius) + "m", {'x [m]' 'y [m]'}, circle.size);
            end
            
            for i = 1:length(obj.speed)
                gamma = 1/(1-(obj.speed(i)^2/obj.c^2)); % Lorentz factor
                % Angle between inner (pointing a middle point of a
                % rod) and outer (pointing ends of a rod) circle radii
                beta = 360/(2*obj.noOfRods*gamma);
                a = obj.radius/cosd(beta); % Outer circle radius

                % Empty array for storing rod's coordinates
                rod = {};
                % Empty array for storing radii pointing the middle of each rod coordinates
                r = {};
                
                % Calculate end-point coordinates for rods and radii
                for j = 1:noOfRods
                    ang = (j-1)*(360/obj.noOfRods);

                    % Calculating coordinates of each of the end points of the rod
                    x1 = a*sind(ang - beta);
                    y1 = a*cosd(ang - beta);
                    x2 = a*sind(ang + beta);
                    y2 = a*cosd(ang + beta);

                    % Calculating coordinates of middle points of the rod
                    x0 = (x1 + x2)/2;
                    y0 = (y1 + y2)/2;

                    % Generating coordinate array for a rod
                    rod(j,:) = {[x1 x2] [y1 y2]};
                    r(j,:) = {[0 x0] [0 y0]};
                end
                
                % Add plots to a tile
                for j = 1:noOfRods
                    layout = addPlot(layout, i, rod(j,:), '-', 'k', 4, "Tangent");
                    layout = addPlot(layout, i, r(j,:), '-.', 'b', 0.5, "Radius");
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
    end
end
        