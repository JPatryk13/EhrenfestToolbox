%% QUANTUMN
% based on given linear speed of an electron and radius of its circular
% path it calcualtes quantum numbers of electron's allowed quantum states.
% 
%   obj = QuantumN(radius, 'speed', speedValue)
%   obj = QuantumN(radius, 'energy', energyValue)
%   obj = QuantumN(radius, Name, Value)
%       constructor, validates input, calculates initial values of the
%       energy (classical kinetic energy) and quantum number. It utilises
%       while loop to determine the range of *allowed quantum numbers.
%       
%       *Each of quantum numbers must provide a corresponding energy level
%       so that all the levels fit in the classical kinetic (initial)
%       energy.
%
%           Input:
%       'radius':           (required), positive number, radius of the
%                           electron's circular path.
%       'speed':            -1 (default*), positive number smaller than the
%                           speed of light, linear speed of an electron in
%                           m/s.
%       'energy':           -1 (default*), positive number, energy of an
%                           electron in joules.
%       'relCorrection':    false (default), boolean; if true, applies
%                           relativistic correction from perturbation
%                           theory to the kinetic energy (quantum number
%                           equation changes as well as it is derived from
%                           the kinetic energy).
%           Output:
%       'obj':              object of the class.
%
%   quantNumbers = getTheList(obj)
%   quantNumbers = getTheList(obj, directQuantNFlag)
%       return list of quantum numbers
%
%           Input:
%       'obj':              object of the class
%       'directQuantNFlag': false (default), if specified by user otherwise
%                           then the function will return row cell array
%                           with the list of quantum numbers as the first
%                           element and the list of quantum numbers
%                           directly calculated from the energy (initial or
%                           fraction of the initial).
%           Output:
%       'quantNumbers':     list of quantum numbers (and optionally direct
%                           quantum numbers list).
%
%   Limitations:
%       N/A
%
%   Examples:
%       Functions quantumNEx1 and quantumNEx2 from Examples.m
%
%   Updates:
%       10/03/2020: Added perturbation theory correction for energy as well
%           as possibility of feeding the function with either speed or
%           kinetic energy. Therefore, input of the constructor was changed
%           - before:(linearSpeed, radius).
%       03/05/2020: Added input parser.
%
%   Use:
%       Such a list (each element) can be fed into WAVEFUNCTION object 
%       to generate allowed wavefunctions of an electron (taking a circular
%       path in a mangetic field).
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, CIRCLE,
%       ENERGYAPPROXIMATION, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE,
%       FINDLIMITS
%
%   Patryk Jesionka, Maciej Makuch, 2019
%%

classdef QuantumN
    properties
        hbar = 1.05*10.^(-34); % Modified Planck's constant
        me = 9.1094*10.^(-31); % Electron rest mass
        c = 3*10^8; % Speed of light in vacuum
        
        % Quantum number value calculated directly from the energy
        directQuantN = []
        % List storing quantum numbers
        quantNumbers = []
    end
    methods
        function obj = QuantumN(radius, varargin)
            % Define default values
            defaultSpeed = -1;
            defaultEnergy = -1;
            defaultRelCorrection = false;
            
            % Validation functions
            validRadius = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x);
            validSpeed = @(x) ((gt(x, 0) && gt(obj.c, x)) || eq(x, -1)) && isreal(x) && isscalar(x) && isnumeric(x);
            validEnergy = @(x) (gt(x, 0) || eq(obj.c, x)) && isreal(x) && isscalar(x) && isnumeric(x) && isfinite(x);
            validRelCorrection = @(x) islogical(x) && isscalar(x);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addRequired(p, 'radius', validRadius);
            addParameter(p, 'speed', defaultSpeed, validSpeed);
            addParameter(p, 'energy', defaultEnergy, validEnergy);
            addParameter(p, 'relCorrection', defaultRelCorrection, validRelCorrection);
            
            parse(p, radius, varargin{:});
            
            % Extract variables from the parser
            radius = p.Results.radius;
            speed = p.Results.speed;
            energy = p.Results.energy;
            relCorrection = p.Results.relCorrection;
            
            % Define energy due to mass
            restMassEnergy = obj.me*obj.c^2;
            
            % Determine whether speed or energy was specified and return an
            % error when both or none were
            if eq(speed, -1) && eq(energy, -1)
                error("Value of speed or energy must be specified.");
            elseif  ~eq(speed, -1) && ~eq(energy, -1)
                error("Only one of the energy and speed values can be specified.");
            elseif eq(speed, -1) 
                % Input is the kinetic energy
                initialEnergy = energy;
            else
                % Input is the speed
                if relCorrection
                    % Considering relativistic approximation
                    initialEnergy = restMassEnergy/sqrt(1-(speed/obj.c)^2) - restMassEnergy;
                else
                    % Considering classical system - no correction
                    initialEnergy = (obj.me*speed^2)/2;
                end
            end
            
            % Energy for a particle on a ring - no perturbation applied
            unperturbedEnergy = @(n) ((n^2)*(obj.hbar^2))/(2*obj.me*(radius^2));
            
            if relCorrection 
                % Relativistically corrected n number
                quantumNFunc = @(E) sqrt((restMassEnergy - restMassEnergy*sqrt(1-2*(E/restMassEnergy)))/unperturbedEnergy(1));
                
                % Correction for relativistic system was from perturbation
                % theory
                correction = @(n) (1/(2*restMassEnergy))*(unperturbedEnergy(n))^2;
            else
                % Quantum number based on the classical energy and circular
                % trajectory radius
                quantumNFunc = @(E) (radius/obj.hbar)*sqrt(2*obj.me*E);
                
                % No correction for classical system
                correction = @(n) 0;
            end
            
            % Obtain quantum number directly from the classical energy
            quantumN = quantumNFunc(initialEnergy);
            
            % Loop through the process as long as quantum number is greater
            % or equal 2
            while quantumN >= 2
                % For high energy input quantumN is complex valued, though,
                % required is  a proper measure giving complex numbers out.
                % Therefore, the loop is carried out with absolute value of
                % the number while the direct quantum number is saved
                % separately.
                obj.directQuantN = [obj.directQuantN quantumN];
                if ~isreal(quantumN)
                    quantumN = abs(quantumN);
                end
                
                if mod(floor(quantumN), 2) == 0
                    % Overwrite the quantum number when its floor is an 
                    % even number
                    quantumN = floor(quantumN);
                elseif floor(quantumN) > 2
                    % Overwrite the quantum number with lower integer when 
                    % its floor is greater than 2 
                    quantumN = floor(quantumN) - 1;
                else
                    % If neither both, break the loop
                    break
                end
                
                % Add quantum number to the list
                obj.quantNumbers = [obj.quantNumbers, quantumN];
                
                % Overwrite energy with the quantum energy subtracted from
                % the initial energy (the value from the past iteration or
                % the classical initial energy)
                initialEnergy = initialEnergy - (unperturbedEnergy(quantumN) - correction(quantumN));
                
                quantumN = quantumNFunc(initialEnergy);
            end
        end
        
        function quantNumbers = getTheList(obj, varargin)
            % Define default values
            defaultDirectQuantNFlag = false;
            
            % Validation functions
            validDirectQuantNFlag = @(x) isscalar(x) && islogical(x);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;
            
            % Adding arguments
            addOptional(p, 'directQuantNFlag', defaultDirectQuantNFlag, validDirectQuantNFlag);
            
            parse(p, varargin{:});
            
            % Extract variables from the parser
            directQuantNFlag = p.Results.directQuantNFlag;
            
            if directQuantNFlag
                % Return the list of quantum numbers with the direct
                % quantum number
                quantNumbers = {obj.quantNumbers obj.directQuantN};
            else
                % Return the list of quantum numbers
                quantNumbers = obj.quantNumbers;
            end
        end
    end
end