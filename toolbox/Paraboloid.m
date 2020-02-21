classdef Paraboloid
    properties
        % Vertical coordinate it should has its tip at
        verticalDisplacement {mustBeNumeric} = 0
        % How flat it should be
        flatness {mustBePositive} = 1
        % How wide it should be
        radius {mustBePositive} = 1
        % Whether it should be directed upwards or downwards
        orientation {mustBeMember(orientation, {'+', '-'})} = '+'
        
        % Structure to be returned with coordinates and dimensions of the
        % plot
        paraboloid = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Paraboloid(verticalDisplacement, flatness, radius, orientation)
            % Input validation
            obj.verticalDisplacement = verticalDisplacement;
            obj.flatness = flatness;
            obj.radius = radius;
            obj.orientation = orientation;
            
            % "Returns 2-D grid coordinates based on the coordinates contained in
            % vectors x and y. X is a matrix where each row is a copy of x, and Y
            % is a matrix where each column is a copy of y."
            [r, ang] = meshgrid(0:0.1:radius, 0:pi/20:2*pi);

            % Based on the generated r and ang, creates coordinate matices 
            x = r.*cos(ang);
            y = r.*sin(ang);
            z = (x.^2)/flatness + (y.^2)/flatness + verticalDisplacement;
            
            obj.paraboloid.coordinates = {x y z};
            
            % Generates size of the plot (tile)
            xlim = [-radius radius];
            ylim = [-radius radius];
            zlim = [verticalDisplacement verticalDisplacement+2*radius];
            
            obj.paraboloid.size = [xlim ylim zlim]; 
        end
        
        function parab = getParaboloid(obj)
            parab = obj.paraboloid;
        end
    end
end

