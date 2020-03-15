% CORONA
%
%   Limitations:
%       Cannot save you from corona.
%
%   Examples:
%       N/A
%
%   Updates:
%       N/A
%   Use:
%       Don't.
%
%   See also:
%       PARABOLOID, SPIRAL, WAVEFUNCTION, PLOT, CIRCLE,
%       ENERGYAPPROXIMATION, WAVE, QUANTUMN
%
%   Patryk Jesionka, Maciej Makuch, 2020

format compact

coronaCases = readmatrix('DailyConfirmedCases.xlsx');

days = coronaCases(:,1) - 43861;
days2 = 0:50;
cases = coronaCases(:,2);

curvefit = fit(days, cases, 'exp2');
func = @(x) curvefit.a*exp(curvefit.b*x) + curvefit.c*exp(curvefit.d*x);

plot(days, cases, 'o', days2, func(days2), '-');
title("Number of COVID-19 cases from 31/01/2020 (The 0 Day)");
legend("Number of confirmed cases", "a*exp(b*days)");
xlabel("days (from 31/01/2020)")
ylabel("number of cases")
for i = 1:7
    text(42+i, func(42+i), int2str(func(42+i)));
end
grid on
