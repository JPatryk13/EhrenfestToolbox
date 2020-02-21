% 
%      INSTRUCTION
% Call in the order:
% 1. Spiral - contructor
% 1a. setCentrePoint (optional)
% 2. setData
% 3. create
% 4. plotSpiral
%
%
classdef Spiral
    properties
        % rmax (later r), rmin: 'starting' and 'ending' radii of the spiral
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

            % Array of radii going from rmax to rmin with a step dr
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
        
        function sprl = getSpiral(obj)
            sprl = obj.spiral;
        end
    end
end