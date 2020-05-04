%% CURRENTDENSITY
% Based on a given radius calculates distribution of the current density
% values across the energy range from 0 to the energy corresponding to the
% speed of light in the vacuum for an electron travelling around a circular
% path. Considers lower speed limitation due to the perturbation if the
% correction is applied.
%
%   obj = CurrentDensity()
%       Constructor, validates user input, preallocates variables and 
%       arrays and predefines anonymous function. Loops through the each 
%       value of the speed specified as a fraction of the speed of light in
%       the vacuum. QUANTUMN is called in each iteration, therefore, the 
%       list is supplied to the inner loops - whether the correction is to
%       be applied or not - then using predefined function for kinetic 
%       energy and current density adds up the values for each quantum 
%       number, thus, adding them to appropriate list. Lastly, the function
%       saves properties into the currentDensity structure. Additionally
%       maximum velocity allowed is calculated and saved in the structure.
%
%           Input:
%       'radius':           (required), radius of the circular path of an
%                           electron.
%       'relCorrection':    false (default), logical value. If set to true,
%                           correction from perturbation theory is applied. 
%       'noOfSamples':      100 (default), number of steps the loop should
%                           take; it can be understood as the quality
%                           factor used in the other functions.
%           Output:
%       'obj':              object of the class
%
%   currentDensity = getCurrentDensity(obj)
%       extract generated data.
%
%           Input:
%       'obj':              object of the class
%           Output:
%       'currentDensity':	structure containing 'coordinates', 'size' and
%                           'maxVel' arrays.
%
%   Limitations:
%       N/A
%
%   Examples:
%       Function currentDensityEx in Examples.m
%
%   Use:
%       The data returned may be used to compare influence of the
%       perturbation correction on the current density of an electron
%       moving around circular path.
%
%   See also:
%       PARABOLOID, CIRCLE, WAVEFUNCTION, PLOT, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE,
%       FINDLIMITS, SPIRAL, MAGNETICFLUX
%
%   Patryk Jesionka, Maciej Makuch 2020
%%

classdef CurrentDensity
    properties
        % Constants
        hbar = 1.0546*10^(-34);         % modified Planck cosntant
        e = 1.6022*10^(-19);            % electron charge
        me = 9.1094*10^(-31);           % electron rest mass
        c = 3*10^8;                     % speed of light in vacuum
        restMassEnergy = (9.1094*10.^(-31))*((3*10^8)^2); % energy due to mass
        
        % Structure to be returned with coordinates and dimensions of the plot
        currentDensity = struct('coordinates', [], 'size', [], 'maxVel', [])
    end
    methods
        function obj = CurrentDensity(radius, varargin)
            % Define default values
            defaultRelCorrection = false;
            defaultNoOfSamples = 100;

            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validRelCorrection = @(x) islogical(x) && isscalar(x);
            validNoOfSamples = @(x) gt(x, 0) && isnumeric(x) && isfinite(x) && isscalar(x) && eq(floor(x), x);

            % Input parser
            p = inputParser;
            p.CaseSensitive = true;

            % Adding arguments
            addRequired(p, 'radius', validRadius);
            addParameter(p, 'relCorrection', defaultRelCorrection, validRelCorrection);
            addParameter(p, 'noOfSamples', defaultNoOfSamples, validNoOfSamples);

            parse(p, radius, varargin{:});

            % Extract variables from the parser
            radius = p.Results.radius;
            relCorrection = p.Results.relCorrection;
            noOfSamples = p.Results.noOfSamples;
            
            % Speed relative to the speed of light - (noOfSamples + 1) due
            % to the limit of speed lower than the speed of light, one
            % additional value of speed not to be considered in the loop
            relVel = (1/(noOfSamples + 1))*obj.c;
            
            % Predefine arrays for coordinates
            CurrDens = zeros(1, noOfSamples);
            EKinetic = zeros(1, noOfSamples);
            
            % Current density and kinetic energy functions
            lorenzFactor = @(v) 1./sqrt(1 - (v/obj.c).^2);
            currentDensity = @(n, r) (obj.e.*obj.hbar.*n)./(4.*pi.*(r.^2).*obj.me);
            kineticEnergy = @(n, r) ((n.^2).*(obj.hbar.^2))./(2.*obj.me.*(r.^2));
            kineticEnergyCorrection = @(n, r) (1/(2*obj.restMassEnergy))*(((n^2)*(obj.hbar^2))/(2*obj.me*(r^2)))^2;
            % Where, v - relative velocity
            %        n - quantum number
            %        r - radius
            
            for i = 1:noOfSamples
                quantumN = QuantumN(radius, 'speed', i*relVel,...
                                            'relCorrection', relCorrection);
                quantumN = getTheList(quantumN, true);
                
                % The case when getTheList function returned complex number
                % as a second parameter means that the energy of an
                % electron exceded limit due to the perturbation correction
                list = quantumN{1};
                directQuantN = quantumN{2};
                if ~isreal(directQuantN)
                    obj.currentDensity.maxVel = i*relVel;
                    break
                end
                
                % Begin the pool with zero-values of current desity and
                % kinetic energy
                Current = 0;
                Ek = 0;

                if relCorrection
                    % Loop through all returned quantum numbers list
                    % applying relativistic correction
                    for j = 1:length(list)
                        % Add the current density due to a particular quantum number to the pool
                        Current = Current + lorenzFactor(relVel).*currentDensity(list(j), radius);
                        % Add the relativistic energy due to a particular quantum number to the pool
                        Ek = Ek + kineticEnergy(list(j), radius) - kineticEnergyCorrection(list(j), radius);
                    end
                else
                    % No correction applied
                    for j = 1:length(list)
                        Current = Current + currentDensity(list(j), radius);
                        Ek = Ek + kineticEnergy(list(j), radius);
                    end
                end
                
                % Add current and kinetic energy pools to the list (as
                % an i^th element). Current density in A*(probability
                % density), kinetic energy in eV
                CurrDens(i) = Current;
                EKinetic(i) = Ek/obj.e;
            end
            
            obj.currentDensity.coordinates = {EKinetic, CurrDens};
            obj.currentDensity.size = [findLimits(EKinetic) findLimits(CurrDens)];
        end
        
        function currentDensity = getCurrentDensity(obj)
            currentDensity = obj.currentDensity;
        end
    end
end