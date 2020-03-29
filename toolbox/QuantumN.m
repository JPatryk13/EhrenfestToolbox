% QUANTUMN - based on given linear speed of an electron and radius of its
% circular path it calcualtes quantum numbers of electron's allowed quantum
% states.
%
%   obj = QuantumN(input, type, radius, relCorrection), constructor, validates input,
%       calculates initial values of the energy (classical kinetic energy)
%       and quantum number. It utilises while loop to determine the range
%       of *allowed quantum numbers.
%       *Each of quantum numbers must provide a corresponding energy level
%       so that all the levels fit in the classical kinetic (initial)
%       energy.
%           Input:
%       'input': positive number, linear speed of an electron or its
%       kinetic energy
%       'type': string of characters, can take either 'speed' or 'energy'
%       value
%       'radius': positive number, radius of the electron's circular path
%       'relCorrection': boolean; if true, applies relativistic correction
%       from perturbation theory to the kinetic energy (quantum number
%       equation changes as well as it is derived from the kinetic energy).
%           Output:
%       'obj': object of the class
%
%   quantNumbers = getTheList(obj), return list of quantum numbers
%           Input:
%       'obj': object of the class
%           Output:
%       'quantNumbers': list of quantum numbers
%
%   Limitations:
%       N/A
%
%   Examples:
%       Let an electron travel with speed v=15000m/s with circular path of
%       the radius r=3?m.
%           quantumN = QuantumN(0.000003, 15000);
%           list = getTheList(quantumN);
%       Output: 390    16     6     4     2
%
%   Updates:
%       10/03/2020: Added perturbation theory correction for energy as well
%       as possibility of feeding the function with either speed or kinetic
%       energy. Therefore, input of the constructor was changed - before:
%       (linearSpeed, radius).
%
%   Use:
%       Such a list (each element) can be fed into WAVEFUNCTION object 
%       to generate allowed wavefunctions of an electron (taking a circular
%       path in a mangetic field).
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, CIRCLE,
%       ENERGYAPPROXIMATION, WAVE, GIF
%
%   Patryk Jesionka, Maciej Makuch, 2019

classdef QuantumN
    properties
        hbar = 1.05*10.^(-34); % Modified Planck's constant
        me = 9.1094*10.^(-31); % Electron rest mass
        c = 3*10^8; % Speed of light in vacuum
        restMassEnergy % Energy due to mass
        
        input {mustBePositive} = 1
        type {mustBeMember(type, {'speed', 'energy'})} = 'energy'
        radius {mustBePositive} = 1
        relCorrection
        
        initialEnergy
        
        quantumN
        quantNumbers = []
    end
    methods
        function obj = QuantumN(input, type, radius, relCorrection)
            % Input validation
            obj.input = input;
            obj.type = type;
            obj.radius = radius;
            if ~islogical(relCorrection)
                error("relCorrection must be logical value!")
            else
                obj.relCorrection = relCorrection;
            end
            
            % Define energy due to mass
            obj.restMassEnergy = obj.me*obj.c^2;
            
            if length(type) == 6 % input is a kinetic energy
                obj.initialEnergy = obj.input;
            else
                % obj.input is the speed value
                if obj.relCorrection % Considering relativistic approximation
                    obj.initialEnergy = obj.restMassEnergy/sqrt(1-(obj.input/obj.c)^2) - obj.restMassEnergy;
                else % Considering classical system - no correction
                    obj.initialEnergy = (obj.me*obj.input^2)/2;
                end
            end
            
            unperturbedEnergy = @(n) ((n^2)*(obj.hbar^2))/(2*obj.me*(obj.radius^2));

            if obj.relCorrection 
                % Relativistically corrected n number
                quantumN = @(E) sqrt((obj.restMassEnergy - obj.restMassEnergy*sqrt(1-2*(E/obj.restMassEnergy)))/unperturbedEnergy(1)); 
                % Correction for relativistic system was from perturbation
                % theory
                correction = @(n) (1/(2*obj.restMassEnergy))*(unperturbedEnergy(n))^2;
            else
                % Quantum number based on the classical energy and circular
                % trajectory radius
                quantumN = @(E) (obj.radius/obj.hbar)*sqrt(2*obj.me*E);
                % Correction for classical system
                correction = @(n) 0;
            end
            
            obj.quantumN = quantumN(obj.initialEnergy);
            
            % Loop through as long as quantum number is greater or equal 2
            while obj.quantumN >= 2
                if mod(floor(obj.quantumN), 2) == 0
                    % Overwrite the quantum number when its floor is an 
                    % even number
                    obj.quantumN = floor(obj.quantumN);
                elseif floor(obj.quantumN) >= 2
                    % Overwrite the quantum number with lower integer when 
                    % its floor is equal or greater than 2 ??? 
                    obj.quantumN = floor(obj.quantumN) - 1;
                else
                    % If neither both, break the loop
                    break
                end
                
                % Add quantum number to the list
                obj.quantNumbers = [obj.quantNumbers, obj.quantumN];
                
                % Overwrite energy with the quantum energy 
                % subtracted from the initial energy (the value from the  
                % past iteration or the classical initial energy)
                obj.initialEnergy = obj.initialEnergy - (unperturbedEnergy(obj.quantumN) - correction(obj.quantumN));
                
                obj.quantumN = quantumN(obj.initialEnergy);
            end
        end
        
        function quantNumbers = getTheList(obj)
            % Return the list of quantum numbers
            quantNumbers = obj.quantNumbers;
        end
    end
end