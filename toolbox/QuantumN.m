% QUANTUMN - based on given linear speed of an electron and radius of its
% circular path it calcualtes quantum numbers of electron's allowed quantum
% states.
%
%   obj = QuantumN(linearSpeed, radius), constructor, validates input,
%       calculates initial values of the energy (classical kinetic energy)
%       and quantum number. It utilises while loop to determine the range
%       of *allowed quantum numbers.
%       *Each of quantum numbers must provide a corresponding energy level
%       so that all the levels fit in the classical kinetic (initial)
%       energy.
%           Input:
%       'linearSpeed': positive number, linear speed of an electron
%       'radius': positive number, radius of the electron's circular path
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
%   Use:
%       Such a list (each element) can be fed into WAVEFUNCTION object 
%       to generate allowed wavefunctions of an electron (taking a circular
%       path in a mangetic field).
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, CIRCLE,
%       ENERGYAPPROXIMATION, WAVE
%
%   Patryk Jesionka, Maciej Makuch, 2019

classdef QuantumN
    properties
        hbar = 1.05*10.^(-34); % Modified Planck's constant
        me = 9.1094*10.^(-31); % Electron rest mass
        linearSpeed {mustBePositive} = 1
        radius {mustBePositive} = 1
        initialEnergy
        quantumN
        quantNumbers = []
    end
    methods
        function obj = QuantumN(linearSpeed, radius)
            % Input validation
            obj.linearSpeed = linearSpeed;
            obj.radius = radius;
            
            % Classical kinetic (initial) energy based on the linear speed
            % of an electron
            obj.initialEnergy = (obj.me*obj.linearSpeed^2)/2;
            % Quantum number based on the classical energy and circular
            % trajectory radius
            obj.quantumN = (obj.radius/obj.hbar)*sqrt(2*obj.me*obj.initialEnergy);
            
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
                obj.initialEnergy = obj.initialEnergy - ((obj.quantumN^2)*(obj.hbar^2))/(2*obj.me*(obj.radius^2));
                
                % Based on the energy value, calculate new quantum number
                obj.quantumN = (obj.radius/obj.hbar)*sqrt(2*obj.me*obj.initialEnergy);
            end
        end
        
        function quantNumbers = getTheList(obj)
            % Return the list of quantum numbers
            quantNumbers = obj.quantNumbers;
        end
    end
end