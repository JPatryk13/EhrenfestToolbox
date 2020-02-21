% WAVEFUNCTION - generating data required to plot a flower-shaped graph of
% the wavefunction. It draws the amplitude in the same plane as the
% circular trajectory of a particle.
%
%   obj = Wavefunction(radius, quantumN), constructor,
%       validates user input, defines base properties of the wavefunction
%       and its domains.
%           Input:
%       'radius': nonnegative number, radius of a circle to create
%       'quantumN': nonnegative, even integer value, quantum number for
%       frequency, period and wavedomain of the function.
%           Output:
%       'obj': object of the class.
%
%   wavefuncion = getWavefunc(obj, time, arithmeticType), validates input,
%       generates wavefunction of the type specified and returns the
%       structure.
%           Input:
%       'obj': object of the class.
%       'time': nonnegative value, wavefunction is time dependent therefore
%       the output coordinates vary with time parameter.
%       'arithmeticType': 'sin' or 'cos' string, specifies which of the
%       component of the wavefunction must be returned with the structure.
%           Output:
%       'wavefunction': structure containing 'coordinates' and 'size'
%       arrays.
%
%   Limitations:
%       Data passed to the functions must be physically sensible - e.g.
%       time cannot be negative, quantum number must be positive, even
%       number.
%
%   Examples:
%       Plot 'sin' component of a wavefunction for the time, t=0 and
%       quantum state described by n=6 for an electron travelling around
%       the circle of radius, r=0.1
%           wavefunctionHandle = Wavefunction(0.1, 6);
%           wavefunction = getWavefunc(wavefunctionHandle, 0, 'sin');
%           plot3(wavefunction.coordinates{1}, wavefunction.coordinates{2},
%                 wavefunction.coordinates{3});
%           axis(wavefunction.size)
%
%   Use:
%       Such a structure ('wavefunction') can be fed into standard MATLAB
%       functions (e.g. plot(), plot3()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, PLOT, QUANTUMN, ENERGYAPPROXIMATION
%
%   Patryk Jesionka, 2019

classdef Wavefunction
    properties
        radius {mustBePositive} % radius of the electron's path
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
            % Input validation
            obj.radius = radius;
            obj.quantumN = quantumN;
            % verification whether quantumN is even or not
            if ~(mod(obj.quantumN, 2) == 0)
                error("Quantum number n must be even!");
            end
            
            % Calcualating the ratio
            obj.rat = obj.hbar/obj.me;
            
            % defining wave properties
            obj.freq = ((obj.quantumN.^2).*obj.rat)./(2*(obj.radius.^2));
            obj.period = (2*pi)/obj.freq;
            obj.amp = (pi*obj.radius)^(-0.5);
            
            % defining domains
            obj.circleDomain = 0:(obj.q*2*pi):(2*pi);
            obj.waveDomain = 0:(obj.q*obj.quantumN*2*pi):(obj.quantumN*2*pi);
        end
        
        function wavefuncion = getWavefunc(obj, time, arithmeticType)
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
            
            % Defining the size of the plot
            pltBounds = obj.radius + obj.amp;
            obj.wavefunc.size = [-pltBounds pltBounds -pltBounds pltBounds -pltBounds pltBounds];
            
            % Returning the structure
            wavefuncion = obj.wavefunc;
        end
    end
end