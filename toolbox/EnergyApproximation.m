classdef EnergyApproximation
    properties
        m = 9.109*10^(-31); % Electron mass
        c = 2.998*10^(8); % Speed of light
        
        min_vel {mustBePositive} = 0.00; % Lower boundary of velocity (relative to the speed o light)
        max_vel {mustBePositive} = 0.99; % Higher boundary of velocity
        step {mustBePositive} = 2; % Velocity resolution relative to 10^5
        
        classical_p % Classical momentum
        relativistic_p % Relativistic momentum
        
        v % Velocity range
        rel_speed % Velocity relative to the speed of light
        exact_value % Exact value of the relativistic kinetic energy
    end
    methods (Access = private)
        function approx = powSeriesApprox(p)
            % 2nd, 3rd and 4th order power series approximations for kinetic energy
            approx_2nd = (obj.m*(obj.c^2)).*(0.5.*((p./(obj.m*obj.c)).^2) - 0.125.*((p./(obj.m*obj.c)).^4));
            approx_3rd = (approx_2nd + (1/16).*((p./(obj.m*obj.c)).^6));
            approx_4th = (approx_3rd - (5/144).*((p./(obj.m*obj.c)).^8));
            
            approx = [approx_2nd approx_3rd approx_4th];
        end
        
        function dev = standardDeviation(p)
            approx = powerSeriesApprox(p);
            
            % Standard deviation of the power series approximations from energy
            dev_2nd_comp1 = (((obj.exact_value.^2)-(approx{1}.^2)).^0.5);
            dev_3rd_comp1 = (((obj.exact_value.^2)-(approx{2}.^2)).^0.5);
            dev_4th_comp1 = (((obj.exact_value.^2)-(approx{3}.^2)).^0.5);
            
            % Standard deviation of energy from the power series approximations
            dev_2nd_comp2 = ((-(obj.exact_value.^2)+(approx{1}.^2)).^0.5);
            dev_3rd_comp2 = ((-(obj.exact_value.^2)+(approx{2}.^2)).^0.5);
            dev_4th_comp2 = ((-(obj.exact_value.^2)+(approx{3}.^2)).^0.5);
            
            % Summation of standard deviation components
            dev_2nd = dev_2nd_comp1 + dev_2nd_comp2;
            dev_3rd = dev_3rd_comp1 + dev_3rd_comp2;
            dev_4th = dev_4th_comp1 + dev_4th_comp2;
            
            dev = [dev_2nd dev_3rd dev_4th];
        end
    end
    methods (Access = public)
        function obj = EnergyApproximation(min_vel, step, max_vel)
            % User input verification
            obj.min_vel = min_vel;
            obj.step = step;
            obj.max_vel = max_vel;
            
            % Calculation - momentum
            obj.relativistic_p = (m*v)./((1-((v./c).^2)).^(0.5)); 
            obj.classical_p = m.*v;

            % Calculation - velocity
            obj.v = (min_vel*c):(step*10^(5)):(max_range*c);
            obj.rel_speed = v./c;
            
            % Calculation - energy
            obj.exact_value = obj.m*(c^2).*(((1+((obj.relativistic_p./(obj.m*c)).^2)).^(0.5))-1);
            
            relativistic_approx = obj.powSeriesApprox(obj.relativistic_p);
            classical_approx = obj.powSeriesApprox(obj.classical_p);
            
            relativistic_dev = standardDeviation(obj.relativistic_p);
            classical_dev = standardDeviation(obj.classical_p);
            
            relativistic_dev_rel = relativistic_dev/exact_value;
            classical_dev_rel = classical_dev/exact_value;
        end
        
        function plotTile(obj)
            % use plot toolbox
        end
    end
end