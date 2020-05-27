%% WAVE
% Generating data required to plot a graph of the wavefunction with
% amplitude in z direction.
%
%   obj = Wave(radius, quantumN)
%   obj = Wave(radius, quantumN, q)
%       Constructor, validates input using CIRCLE and WAVEFUCNTION classes,
%       generates circle coordinates (x and y of use) and initialises
%       wavefunction object.
%
%           Input:
%       'radius':       (required), nonnegative number, radius of the
%                       circle the wavefunction will be plotted on.
%       'quantumN':     (required), nonnegative, even integer value,
%                       quantum number for frequency, period and wavedomain
%                       of the function.
%       'q':            360 (default), must be positive integer, quality
%                       factor (see Updates).
%           Output:
%       'obj':          object of the class.
%
%   wave = getWave(obj, time)
%   wave = getWave(obj, time, arithmeticType)
%       Validates input using WAVEFUNCTION class, assigns data to the
%       structure.
%
%           Input:
%       'obj':          object of the class.
%       'time':         0 (default), nonnegative value, wavefunction is
%                       time dependent therefore the output coordinates
%                       vary with time parameter.
%       'arithmeticType': 'sin' (default), 'sin' or 'cos' string, specifies
%                       which of the component of the wavefunction must be
%                       returned with the structure.
%           Output:
%       'wave':         structure containing 'coordinates' and 'size'
%                       arrays.
%
%   Limitations:
%       Such a wavefunction has it's centre point at [0 0 0] and currently
%       it cannot be modified.
%
%   Examples:
%       Function waveEx from Examples.m
%
%   Updates:
%       29/04/2020: Added input parser.
%
%   Use:
%       Such a structure ('wave') can be fed into standard MATLAB functions
%       (e.g. plot3()). However, its purpose is to input data into
%       the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, PLOT, QUANTUMN, ENERGYAPPROXIMATION,
%       WAVEFUNCTION, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE, FINDLIMITS,
%       CURRENTDENSITY, MAGNETICFLUX
%
%   Patryk Jesionka, 2020
%%

classdef Wave
    properties
        wavefunctionHandle      % Handle for the Wavefunction object
        
        % Structure containing coordinates and dimensions of the circle the
        % wavefunction should be plotted on
        circle = struct('coordinates', [], 'size', [])
        
        % Structure to be returned with coordinates and dimensions of
        % the plot
        wave = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Wave(radius, quantumN, varargin)
            % Define default values
            defaultQ = 360;
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'radius');
            addRequired(p, 'quantumN');
            addOptional(p, 'q', defaultQ);
            
            parse(p, radius, quantumN, varargin{:});
            
            % Extract variables from the parser
            radius = p.Results.radius;
            quantumN = p.Results.quantumN;
            q = p.Results.q;
            
            % Use CIRCLE and WAFUNCTION class objects to verify input and
            % create circle as a base for the wavefunction and to
            % initialise wavefunction object
            circleHandle = Circle(radius, 'q', q);
            obj.circle = getCircle(circleHandle);
            obj.wavefunctionHandle = Wavefunction(radius, quantumN, 'q', q);
            
            % Set x and y axis limits of the plot
            xlim = findLimits(obj.circle.coordinates{1});
            ylim = findLimits(obj.circle.coordinates{2});
            obj.wave.size(1:4) = [xlim ylim];
        end
        
        function wave = getWave(obj, varargin)
            % Define default values
            defaultTime = 0;
            defaultArithmeticType = 'sin';
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addParameter(p, 'time', defaultTime);
            addParameter(p, 'arithmeticType', defaultArithmeticType);
            
            parse(p, varargin{:});
            
            % Extract variables from the parser
            time = p.Results.time;
            arithmeticType = p.Results.arithmeticType;
            
            % Return coordinates and the size from the wavefunction object
            wavefunction = getWavefunc(obj.wavefunctionHandle,...
                'time', time,...
                'arithmeticType', arithmeticType,...
                'amplitudeAxes', 'z');
            
            % Put data in the structure
            obj.wave.coordinates = {obj.circle.coordinates{1} obj.circle.coordinates{2} wavefunction.coordinates{3}};
            obj.wave.size(5:6) = wavefunction.size(5:6);
            
            % Return the data
            wave = obj.wave;
        end
    end
end

