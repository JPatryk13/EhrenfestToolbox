classdef Wavefunction
    properties
        radius {mustBeNumeric} % radius of the electron's path
        quantumN {mustBeInteger} % quantum number n
        hbar = 1.05*10.^(-34); % modified Planck's constant
        me = 9.1094*10.^(-31); % electron rest mass
        rat % ratio of Planck's constant and electron rest mass
        freq % frequency of the wave
        q = 0.005; % quality factor, the lower it is the more steps plotter takes (time vs quality)
        circleDomain % range from 0 to 2*pi (generating a circle)
        waveDomain % range from 0 to 2*n*pi (generating a weve)
        amp % amplitude of the wave (normalisation constant)
        period % time of the one escillation
        % Structure to be returned with coordinates and dimensions of the plot
        wavefunc = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Wavefunction(radius, quantumN)
            obj.radius = radius;
            obj.quantumN = quantumN;
            % verification whether n_ (quantum no.) is even or not
            if ~(mod(obj.quantumN, 2) == 0)
                error("Quantum number n must be even!");
            end
            
            obj.rat = obj.hbar/obj.me;
            
            % defining wave properties
            obj.freq = ((obj.quantumN.^2).*obj.rat)./(2*(obj.radius.^2));
            obj.period = (2*pi)/obj.freq;
            obj.amp = (pi*obj.radius)^(-0.5);
            
            % defining domains
            obj.circleDomain = 0:(obj.q*2*pi):(2*pi);
            obj.waveDomain = 0:(obj.q*obj.quantumN*2*pi):(obj.quantumN*2*pi);
        end
        
        function wavefunc = plotWavefunc(obj, time, arithmeticType)
            % Input validation
            if time < 0
                error("Time must be a positive value!");
            elseif ~ismember(arithmeticType, {'sin', 'cos'})
                error("arithmeticType must be either 'sin' or 'cos'!");
            end
            
            % Generating coordinates
            if arithmeticType == 'sin'
                xSin = obj.radius*cos(obj.circleDomain) + obj.amp.*sin(obj.waveDomain).*cos(obj.circleDomain).*cos(obj.freq.*time);
                ySin = obj.radius*sin(obj.circleDomain) + obj.amp.*sin(obj.waveDomain).*sin(obj.circleDomain).*cos(obj.freq.*time);
                zSin = transpose(-obj.amp.*sin(obj.freq.*time).*ones(length(xSin), 1));
                obj.wavefunc.coordinates = {xSin ySin zSin};
            else
                xCos = obj.radius*cos(obj.circleDomain) + obj.amp.*cos(obj.waveDomain).*cos(obj.circleDomain).*cos(obj.freq.*time);
                yCos = obj.radius*sin(obj.circleDomain) + obj.amp.*cos(obj.waveDomain).*sin(obj.circleDomain).*cos(obj.freq.*time);
                zCos = transpose(-obj.amp.*sin(obj.freq.*time).*ones(length(xCos), 1));
                obj.wavefunc.coordinates = {xCos yCos zCos};
            end
            
            % Defining size of the plot
            pltBounds = obj.radius + obj.amp;
            obj.wavefunc.size = [-pltBounds pltBounds -pltBounds pltBounds -pltBounds pltBounds];
            
            % Returning structure
            wavefunc = obj.wavefunc;
        end
    end
end