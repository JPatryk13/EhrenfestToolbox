% WAVE - generating data required to plot a graph of the wavefunction with
% amplitude in z direction.
%
%   obj = Wave(radius, quantumN, q), constructor, validates input using CIRCLE
%       and WAVEFUCNTION classes, generates circle coordinates (x and y of
%       use) and initialises wavefunction object.
%           Input:
%       'radius': nonnegative number, radius of the circle the wavefunction
%       will be plotted on.
%       'quantumN': nonnegative, even integer value, quantum number for
%       frequency, period and wavedomain of the function.
%       'q': must be positive integer, quality factor (see Updates).
%           Output:
%       'obj': object of the class.
%
%   wave = getWave(obj, time, arithmeticType), validates input using
%       WAVEFUNCTION class, assigns data to the structure.
%           Input:
%       'obj': object of the class.
%       'time': nonnegative value, wavefunction is time dependent therefore
%       the output coordinates vary with time parameter.
%       'arithmeticType': 'sin' or 'cos' string, specifies which of the
%       component of the wavefunction must be returned with the structure.
%           Output:
%       'wave': structure containing 'coordinates' and 'size' arrays.
%
%   Limitations:
%       Such a wavefunction has it's centre point at [0 0 0] and currently
%       in cannot be modified.
%
%   Examples:
%       Plot a wavefunction of an electron which energy is described by
%       n=8' It is travelling around the circular path of radius, r=0.1m.
%       Desired type of the wavefunction is cosine at time, t=0.
%           waveHandle = Wave(0.1, 8);
%           wave = getWave(waveHandle, 0, 'cos');
%           plot3(wave.coordinates{1}, wave.coordinates{2},
%                 wave.coordinates{3});
%           axis(wave.size)
%
%   Use:
%       Such a structure ('wave') can be fed into standard MATLAB functions
%       (e.g. plot3()). However, its purpose is to input data into
%       the PlotToolbox's PLOT.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, PLOT, QUANTUMN, ENERGYAPPROXIMATION,
%       WAVEFUNCTION
%
%   Patryk Jesionka, 2020

classdef Wave
    properties
        wavefunctionHandle % Handle for the Wavefunction object
        % Structure containing coordinates and dimensions of the circle the
        % wavefunction should be plotted on
        circle = struct('coordinates', [], 'size', [])
        % Structure to be returned with coordinates and dimensions of the plot
        wave = struct('coordinates', [], 'size', [])
    end
    methods
        function obj = Wave(radius, quantumN, q)
            % Usse CIRCLE and WAFUNCTION class objects to verify input and
            % create circle as a base for the wavefunction and to
            % initialise wavefunction object
            circleHandle = Circle(radius, [0 0 0], 'z', q);
            obj.circle = getCircle(circleHandle);
            obj.wavefunctionHandle = Wavefunction(radius, quantumN, q);
            
            % Set x and y axis limits of the plot
            obj.wave.size(1:4) = [-radius radius -radius radius].*1.5;
        end
        
        function wave = getWave(obj, time, arithmeticType)
            % Return coordinates and the size from the wavefunction object
            wavefunction = getWavefunc(obj.wavefunctionHandle, time, arithmeticType, 'z');
            
            % Put data in the structure
            obj.wave.coordinates = {obj.circle.coordinates{1} obj.circle.coordinates{2} wavefunction.coordinates{3}};
            obj.wave.size(5:6) = wavefunction.size(5:6);
            
            % Return the data
            wave = obj.wave;
        end
    end
end

