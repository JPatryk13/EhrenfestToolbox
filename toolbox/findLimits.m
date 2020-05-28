%% FINDLIMITS
% Finding maximum and minimum values in a given array.
%
%   lim = findLimits(coordinate)
%   lim = findLimits(coordinate, displacement)
%
%           Input:
%       'coordinate':       (required), array with numeric values.
%       'displacement':     0 (default), shift that should be applied to
%                           the 'coordinate' array.
%           Output:
%       'lim':              array with limits, [min max].
%
%   Limitations:
%       Does not support multiple array-displacement pair input.
%
%   Examples:
%       A = [-1 0 1 2 3 4 5 6];
%       dA = -2;
%       findLimits(A, dA)
%
%   Use:
%       Extracting min-max limits pair from a numeric array.
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, PLOT, GIF, DAVIDOVICRODS,
%       CHANGENOTATIONTYPE, CURRENTDENSITY, MAGNETICFLUX
%
%   Patryk Jesionka, 2020
%%

function lim = findLimits(coordinate, varargin)
    % Define default values
    defaultDisplacement = 0;
            
    % Validation functions
    validCoordinate = @(x) all(isnumeric(x)) && all(all(isfinite(x)));
    validDisplacement = @(x) isnumeric(x) && isfinite(x) && isreal(x) && isscalar(x);
            
    % Input parser
    p = inputParser;
    p.CaseSensitive = true;
            
    % Adding arguments
    addRequired(p, 'coordinate', validCoordinate);
    addOptional(p, 'displacement', defaultDisplacement, validDisplacement);
            
    parse(p, coordinate, varargin{:});
            
    % Extract variables from the parser
    a = p.Results.coordinate;
    da = p.Results.displacement;
    
    % Finding minimum and maximum values in the array; double min/max
    % function to get rid of the multiplicity of limits.
    if min(min(a)) == max(max(a))
        % If min = max then return inf values - they set the axis() to
        % default (optimal) value
        lim = [-inf inf];
    else
        lim = [min(min(a)) + da, max(max(a)) + da];
    end
end