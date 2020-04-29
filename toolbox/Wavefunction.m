% WAVEFUNCTION - generating data required to plot a flower-shaped graph of
% the wavefunction. It draws the amplitude in the same plane as the
% circular trajectory of a particle.
%
%   obj = Wavefunction(radius, quantumN)
%   obj = Wavefunction(radius, quantumN, q)
%       constructor, validates user input, defines base properties of the
%       wavefunction and its domains.
%
%           Input parameters' names:
%       'radius':       (required), nonnegative number, radius of a circle
%                       to create.
%       'quantumN':     (required), nonnegative, even integer value,
%                       quantum number for frequency, period and wavedomain
%                       of the function.
%       'q':            360 (default), must be positive integer, quality 
%                       factor (see Updates).
%           Output:
%       'obj':          object of the class.
%
%   wavefuncion = getWavefunc(obj, time)
%   wavefuncion = getWavefunc(obj, time, Name, Value) 
%       validates input, generates wavefunction of the type specified and 
%       returns the structure.
%
%           Input parameters' names:
%       'obj':              object of the class.
%       'time':             (required), nonnegative value, wavefunction is
%                           time dependendent therefore the output
%                           coordinates vary with time parameter.
%       'arithmeticType':   'sin' (default), 'sin' or 'cos' string,
%                           specifies which of the component of the
%                           wavefunction must be returned with the
%                           structure.
%       'amplitudeAxes':    'xy' (required), which type of the data should
%                           be returned (see Updates)
%           Output:
%       'wavefunction':     structure containing 'coordinates' and 'size'
%                           arrays.
%
%   Limitations:
%       Data passed to the functions must be physically sensible - e.g.
%       time cannot be negative, quantum number must be positive, even
%       number.
%
%   Examples:
%       Function wavefunctionEx in Examples.m
%
%   Updates:
%       27/02/2020: Added new input parameter in getWavefunc(),
%           amplitudeAxes - can take either 'xy' or 'z' value. It decides
%           on the type of the returned data. In case of 'xy' the
%           flower-shaped wavefunction coordinates are returned, otherwise
%           x and y coordinate arrays are ampty but the z coordinate;
%           it contains information on wave domain. It allows to plot
%           wavefunction with aplitude in the z-axis.
%       01/03/2020: Added quality factor (constructor). Need for unifying
%           CIRCLE, WAVEFUNCTION and  SPIRAL classes' output number of
%           steps replaced with freedom of choice.
%       29/04/2020: Added input parser.
%
%   Use:
%       Such a structure ('wavefunction') can be fed into standard MATLAB
%       functions (e.g. plot(), plot3()). However, its purpose is to input
%       data into the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, PLOT, QUANTUMN, ENERGYAPPROXIMATION,
%       WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE, FINDLIMITS
%
%   Patryk Jesionka, 2019

classdef Wavefunction
    properties
        radius {mustBePositive} % Radius of the electron's path
        
        hbar = 1.05*10.^(-34); % Modified Planck's constant
        me = 9.1094*10.^(-31); % Electron rest mass
        rat % Ratio of Planck's constant and electron rest mass
        
        freq % Frequency of the wave
        amp % Amplitude of the wave (normalisation constant)
        period % Time of the one escillation
        
        circleDomain % Range from 0 to 2*pi (generating a circle)
        waveDomain % Range from 0 to 2*n*pi (generating a weve)
        
        % Structure to be returned with coordinates and dimensions of the plot
        wavefunc = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Wavefunction(radius, quantumN, varargin)
            % Define default values
            defaultQ = 360;
            
            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validQuantumN = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && eq(mod(x, 2), 0);
            validQ = @(x) validRadius(x) && eq(x, floor(x));
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'radius', validRadius);
            addRequired(p, 'quantumN', validQuantumN);
            addOptional(p, 'q', defaultQ, validQ);
            
            parse(p, radius, quantumN, varargin{:});
            
            % Extract variables from the parser
            obj.radius = p.Results.radius;
            quantumN = p.Results.quantumN;
            q = p.Results.q;
                        
            % Calcualating the ratio
            obj.rat = obj.hbar/obj.me;
            
            % defining wave properties
            obj.freq = ((quantumN.^2).*obj.rat)./(2*(obj.radius.^2));
            obj.period = (2*pi)/obj.freq;
            obj.amp = (pi*obj.radius)^(-0.5);
            
            % defining domains
            obj.circleDomain = 0:(2*pi/q):(2*pi);
            obj.waveDomain = 0:(quantumN*2*pi/q):(quantumN*2*pi);
        end
        
        function wavefuncion = getWavefunc(obj, time, varargin)
            % Define default values
            defaultArithmeticType = 'sin';
            defaultAmplitudeAxes = 'xy';
            
            % Validation functions
            validTime = @(x) ge(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validArithmeticType = @(x) ischar(x) && ismember(x, {'sin', 'cos'});
            validAmplitudeAxes = @(x) ischar(x) && ismember(x, {'xy', 'z'});
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'time', validTime);
            addParameter(p, 'arithmeticType', defaultArithmeticType, validArithmeticType);
            addParameter(p, 'amplitudeAxes', defaultAmplitudeAxes, validAmplitudeAxes);
            
            parse(p, time, varargin{:});
            
            % Extract variables from the parser
            time = p.Results.time;
            arithmeticType = p.Results.arithmeticType;
            amplitudeAxes = p.Results.amplitudeAxes;
            
            % Generating coordinates. Return coordinates of the
            % flower-shaped wavefunction if the 'amplitudeAxes' is
            % specified as 'xy'; otherwise return empty arrays at x and y
            % coordinates and wavefunction with amplitude in z direction.
            if eq(amplitudeAxes, 'xy')
                if eq(arithmeticType, 'sin')
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
                
                % Sets optimal size of the plot
                lim = findLimits(obj.wavefunc.coordinates{1});
            else
                if eq(arithmeticType, 'sin')
                    xSin = [];
                    ySin = [];
                    zSin = obj.amp.*sin(obj.waveDomain).*cos(obj.freq.*time);
                    obj.wavefunc.coordinates = {xSin ySin zSin};
                else
                    xCos = [];
                    yCos = [];
                    zCos = obj.amp.*cos(obj.waveDomain).*cos(obj.freq.*time);
                    obj.wavefunc.coordinates = {xCos yCos zCos};
                end
                % Sets optimal size of the plot
                lim = findLimits(obj.wavefunc.coordinates{3});
            end
            
            obj.wavefunc.size = [lim lim lim];
            
            % Returning the structure
            wavefuncion = obj.wavefunc;
        end
    end
end