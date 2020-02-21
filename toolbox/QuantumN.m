classdef QuantumN
    properties
        hbar = 1.05*10.^(-34); % modified Planck's constant
        me = 9.1094*10.^(-31); % electron rest mass
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
            
            obj.initialEnergy = (obj.me*obj.linearSpeed^2)/2; % classical energy
            obj.quantumN = (obj.radius/obj.hbar)*sqrt(2*obj.me*obj.initialEnergy);
            
            while obj.quantumN >= 2
                if mod(floor(obj.quantumN), 2) == 0
                    obj.quantumN = floor(obj.quantumN);
                elseif floor(obj.quantumN) >= 2
                    obj.quantumN = floor(obj.quantumN) - 1;
                else
                    break
                end
                
                obj.quantNumbers = [obj.quantNumbers, obj.quantumN];
                
                obj.initialEnergy = obj.initialEnergy - ((obj.quantumN^2)*(obj.hbar^2))/(2*obj.me*(obj.radius^2));
                
                obj.quantumN = (obj.radius/obj.hbar)*sqrt(2*obj.me*obj.initialEnergy);
            end
        end
        
        function quantNumbers = getTheList(obj)
            quantNumbers = obj.quantNumbers;
        end
    end
end