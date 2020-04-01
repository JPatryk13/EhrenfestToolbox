% CORONA
%
%   Limitations:
%       Cannot save you from corona.
%
%   Use:
%       Don't.
%
%   Patryk Jesionka, Maciej Makuch, 2020

format compact

fileName = 'DailyConfirmedCases.xlsx';

% Remove old data file if present
if exist(fileName, 'file') == 2
  delete(fileName);
end

% Download fresh data from URL
url = ["https://www.arcgis.com/sharing/rest/content/items/e5fd11150d274bebaaf8fe2a7a2bda11/data"];
websave(fileName, url);

% Get the data
coronaCases = readmatrix(fileName);

% Extract data to arrays
days = coronaCases(:,1) - 43861;
cases = coronaCases(:,3);

days1 = 1:57; % day 56: 27/03/2020
func1 = curveFit(40, coronaCases); % 10/03/2020
days2 = 1:(length(cases) - 1);
func2 = curveFit(55, coronaCases); % 25/03/2020

% Plot data and the curve
plot(days, cases, 'kd',...
     days1, func1(days1), 'r-',...
     days2, func2(days2), 'g-');

% Set title
title("Number of COVID-19 cases from 31/01/2020 (The 0 Day)");

% Add legend
lgd = legend("Number of confirmed cases",...
             "Exponential growth - approximation" + newline + "based on the first 40 days",...
             "Approximation based on the first 55" + newline + "days of breakdown");
lgd.Location = 'northwest';

% Add labels
xlabel("days (from 31/01/2020)");
ylabel("number of cases");

% Display the most recent sample values
text(length(days)-1, cases(length(days)),...
     int2str(cases(length(days))) + "  ", 'HorizontalAlignment', 'right');
text(length(days1), func1(length(days1)),...
     "(" + int2str(func1(length(days1))) + ", 27/03/2020) ", 'HorizontalAlignment', 'right');
text(length(days2), func2(length(days2)),...
    int2str(func2(length(days2))) + "  ", 'HorizontalAlignment', 'right');

% Turn the grid on
grid on

function func = curveFit(day, coronaCases)
    % Extract first 'day' samples - for curve fitting
    daysFit = coronaCases(1:day,1) - 43861;
    casesFit = coronaCases(1:day,3);

    % Set the line of best fit
    curvefit = fit(daysFit, casesFit, 'exp2');
    func = @(x) curvefit.a*exp(curvefit.b*x) + curvefit.c*exp(curvefit.d*x);
end
