%% ENERGYAPPROXIMATION
% based on a given speed range and speed resolution calculate: energy exact
% value, energy approximation, approximation deviation from the exact
% value, deviation relative to the exact value. Thus, given is choice of 
% the output of the get-function. 
% 
%   obj = EnergyApproximation(minVel, maxVel)
%   obj = EnergyApproximation(minVel, maxVel, step)
%       constructor, validates user input. Calculates relative velocity -
%       to be treated as horizontal (x) coordinate and exact, approximated
%       values of the energy as well as relative/non-relative deviation
%       from the exact energy value - to be treated as the vertical
%       coordinate.
%
%           Input:
%       'minVel', 'maxVel':     (required), the lowest and the highest
%                               value of the velocity. It is expressed as a
%                               fraction of the speed of light.
%       'step':                 0.1 (default), resolution of the velocity
%                               range expressed as a fraction of the speed
%                               of light. While creating an array of speed
%                               values treated os a spacing between them.
%           Output:
%       'obj':                  object of the class.
%
%   energyApproximation = getEnergyApproximation(obj, energy)
%   energyApproximation = getEnergyApproximation(obj, energy, Name, Value)
%       ...
%
%           Input:
%       'obj':                  object of the class.
%       'energy':               (required), the data to be returned by the
%                               function - 'exact', 'approximation', etc.
%       'model':                'classical' (default)
%       'order':                2 (default)
%                               Both 'model' and 'order' are related to
%                               data returned by the function, e.g. 
%                               model=relativistic, order=2 with 
%                               energy=approximation will return
%                               2nd order relativistic correction the
%                               energy.
%           Output:
%       'energyApproximation':  structure containing 'coordinates' and
%                               'size' arrays.
%
%   approx = powSeriesApprox(obj, p)
%       (private method) approximates energy using power series based on
%       the momentum 'p'.
%
%           Input:
%       'obj':                  object of the class.
%       'p':                    (required), momentum of an electron.
%           Output:
%       'approx':               cell array with 1st, 2nd and 3rd order
%                               approximation of the energy.
%
%   dev = standardDeviation(obj, p)
%       (private method) calculates standard deviation of an approximation 
%       - given by the power series with momentum 'p' input - from the
%       exact energy value.
%
%           Input:
%       'obj':                  object of the class.
%       'p':                    (required), momentum of an electron.
%           Output:
%       'dev':                  cell array with the standard deviation of
%                               the 1st, 2nd and 3rd order approximation.
%
%   Limitations:
%       N/A
%
%   Examples:
%       Function energyApproximationEx from Examples.m
%
%   Updates:
%       04/05/2020: Added input parser.
%
%   Use:
%       Returned data can be used to plot comparison between different
%       order approximation of the energy of an electron treated classicaly
%       or relativistically. The data allow to show how the approximate
%       values deviate from the exact value. Ideally to be used with PLOT
%       class, however, built-in MATLAB classes are not unrecommended.
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, CIRCLE,
%       QUANTUMN, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE, FINDLIMITS
%
%   Patryk Jesionka, Maciej Makuch, 2019
%%

classdef EnergyApproximation
    properties
        % Constants
        m = 9.109*10^(-31); % Electron mass
        c = 2.998*10^(8); % Speed of light
        
        % Velocity relative to the spee of light
        relVel
        
        % Energy exact value
        exactValue

        % Classical and relativistic approximation of the energy value
        classicalApprox
        relativisticApprox

        % Classical and relativistic energy deviation from the exact
        % value
        classicalDev
        relativisticDev

        % Classical and relativistic energy deviation from the exact
        % value relative to the exact value
        classicalDevRel
        relativisticDevRel
        
        energyApproximation = struct('coordinates', [], 'size', [])
    end
    methods (Access = private)
        function approx = powSeriesApprox(obj, p)
            % 2nd, 3rd and 4th order power series approximations for kinetic energy
            approx_2nd = (obj.m*(obj.c^2)).*(0.5.*((p./(obj.m*obj.c)).^2) - 0.125.*((p./(obj.m*obj.c)).^4));
            approx_3rd = (approx_2nd + (1/16).*((p./(obj.m*obj.c)).^6));
            approx_4th = (approx_3rd - (5/144).*((p./(obj.m*obj.c)).^8));
            
            approx = {approx_2nd approx_3rd approx_4th};
        end
        
        function dev = standardDeviation(obj, p)
            approx = powSeriesApprox(obj, p);
            
            % Standard deviation of the power series approximations from energy
            dev_2nd_comp1 = (((obj.exactValue.^2) - (approx{1}.^2)).^0.5);
            dev_3rd_comp1 = (((obj.exactValue.^2) - (approx{2}.^2)).^0.5);
            dev_4th_comp1 = (((obj.exactValue.^2) - (approx{3}.^2)).^0.5);
            
            % Standard deviation of energy from the power series approximations
            dev_2nd_comp2 = ((-(obj.exactValue.^2) + (approx{1}.^2)).^0.5);
            dev_3rd_comp2 = ((-(obj.exactValue.^2) + (approx{2}.^2)).^0.5);
            dev_4th_comp2 = ((-(obj.exactValue.^2) + (approx{3}.^2)).^0.5);
            
            % Summation of standard deviation components
            dev_2nd = dev_2nd_comp1 + dev_2nd_comp2;
            dev_3rd = dev_3rd_comp1 + dev_3rd_comp2;
            dev_4th = dev_4th_comp1 + dev_4th_comp2;
            
            dev = {dev_2nd dev_3rd dev_4th};
        end
    end
    methods (Access = public)
        function obj = EnergyApproximation(minVel, maxVel, varargin)
            % Define default values
            defaultStep = 0.1; % relative to the speed of light
            
            % Validation functions
            validVel = @(x) gt(x, 0) && gt(1, x) && isreal(x) && isscalar(x) && isnumeric(x);
            validStep = @(x) gt(x, 0) && gt(1, x) && isnumeric(x) && isreal(x) && isscalar(x);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'minVel', validVel);
            addRequired(p, 'maxVel', validVel);
            addOptional(p, 'step', defaultStep, validStep);
            
            parse(p, minVel, maxVel, varargin{:});
            
            % Extract variables from the parser
            minVel = p.Results.minVel;
            maxVel = p.Results.maxVel;
            step = p.Results.step;
            
            % Calculation - velocity
            v = (minVel*obj.c):(step*obj.c):(maxVel*obj.c);
            obj.relVel = v./obj.c;
            
            % Calculation - momentum
            relativisticP = (obj.m*v)./((1-((v./obj.c).^2)).^(0.5)); 
            classicalP = obj.m.*v;
            
            % Calculation - energy exact value
            obj.exactValue = obj.m*(obj.c^2).*(((1+((relativisticP./(obj.m*obj.c)).^2)).^(0.5))-1);
            
            % Classical and relativistic approximation of the energy value
            obj.classicalApprox = powSeriesApprox(obj, classicalP);
            obj.relativisticApprox = powSeriesApprox(obj, relativisticP);
            
            % Classical and relativistic energy deviation from the exact
            % value
            obj.classicalDev = standardDeviation(obj, classicalP);
            obj.relativisticDev = standardDeviation(obj, relativisticP);
            
            % Classical and relativistic energy deviation from the exact
            % value relative to the exact value
            obj.classicalDevRel = { obj.classicalDev{1}./obj.exactValue,...
                                    obj.classicalDev{2}./obj.exactValue,...
                                    obj.classicalDev{3}./obj.exactValue};
            obj.relativisticDevRel = { obj.relativisticDev{1}./obj.exactValue,...
                                       obj.relativisticDev{2}./obj.exactValue,...
                                       obj.relativisticDev{3}./obj.exactValue};
        end
        
        function energyApproximation = getEnergyApproximation(obj, energy, varargin)
            % Define default values
            defaultModel = 'classical';
            defaultOrder = 2;
            
            % Validation functions
            validEnergy = @(x) ismember(x, {'exact' 'approximation' 'deviation' 'deviationRel'});
            validModel = @(x) ismember(x, {'classical' 'relativistic'});
            validOrder = @(x) ismember(x, [2 3 4]);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'energy', validEnergy);
            addParameter(p, 'model', defaultModel, validModel);
            addParameter(p, 'order', defaultOrder, validOrder)
            
            parse(p, energy, varargin{:});
            
            % Extract variables from the parser
            energy = p.Results.energy;
            model = p.Results.model;
            order = p.Results.order;
            
            xLim = findLimits(obj.relVel);
            
            % Use 'order' as indices
            order = order - 1;
            
            % Return approprate coordinate-pair and size arrays
            % length(energy) = 5(exact), 13(approximation), 9(deviation),
            %                  12(deviationRel)
            % length(model) = 9(classical), 12(relativistic)
            switch length(energy)
                case 5
                    energyApproximation.coordinates = {obj.relVel obj.exactValue};
                    energyApproximation.size = [xLim findLimits(obj.exactValue)];
                case 13
                    if length(model) == 9
                        energyApproximation.coordinates = {obj.relVel obj.classicalApprox{order}};
                        energyApproximation.size = [xLim findLimits(obj.classicalApprox{order})];
                    else
                        energyApproximation.coordinates = {obj.relVel obj.relativisticApprox{order}};
                        energyApproximation.size = [xLim findLimits(obj.relativisticApprox{order})];
                    end
                case 9
                    if length(model) == 9
                        energyApproximation.coordinates = {obj.relVel obj.classicalDev{order}};
                        energyApproximation.size = [xLim findLimits(obj.classicalDev{order})];
                    else
                        energyApproximation.coordinates = {obj.relVel obj.relativisticDev{order}};
                        energyApproximation.size = [xLim findLimits(obj.relativisticDev{order})];
                    end
                otherwise
                    if length(model) == 9
                        energyApproximation.coordinates = {obj.relVel obj.classicalDevRel{order}};
                        energyApproximation.size = [xLim findLimits(obj.classicalDevRel{order})];
                    else
                        energyApproximation.coordinates = {obj.relVel obj.relativisticDevRel{order}};
                        energyApproximation.size = [xLim findLimits(obj.relativisticDevRel{order})];
                    end
            end
        end
    end
end