% CHANGENOTATIONTYPE - 
%
%   strNumber = changeNotationType(number, type)
%           Input:
%       'number': float value that meant to be converted to string and
%       displayed with scientific notation.
%       'type': char, either 's' or 'e' - decides on whether scientific or 
%       e notation should be used 
%           Output:
%       'strNumber': converted number as a string of characters
%
%   Limitations:
%       Rounds output value to two decimal places. Can handle numbers from
%       the range order from 10^-40 to 10^40.
%
%   Examples:
%       number = 329863741;
%       changeNotationType(number, 'scientific')
%
%   Use:
%       Converting large (e.g. 329863741) or small (e.g. 0.00000000023)
%       into more readable form (3.29*10^8, 3.29e8, 2.3*10^-10, 2.3e-10).
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, PLOT, GIF, DAVIDOVICRODS
%
%   Patryk Jesionka, 2020
%
% Note: haven't tested for surfaces plotting

function strNumber = changeNotationType(number, type)
    % Define prefix based on the type value
    if type == 's'
        preffix = "\times10^";
    elseif type == 'e'
        preffix = "e";
    else % User input validaiton
        error("type should be a character value, 's' or 'e'.");
    end
    if ~isnumeric(number)
        error("number must be numberic.");
    end
    if (number < 10 && number >= 0.1) || (number > -10 && number <= -0.1) || (number == 0)
        strNumber = string(number);
    else
        for n = -40:40
            if abs(number) >= 10^n
                strNumber = string(number/(10^n)) + preffix + string(n);
            end
        end
    end
end

