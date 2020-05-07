format compact

% davidovicStarDisk

% diskDef

% wavefunction

% correctionOrders
% correctionOrdersLog
% correctionDeviationLog
% correctionDeviationRel2ndOrder
% classVsRelApprox

radius = 35*10^(-9);        % 35nm
energy = 10*1.602*10^-19;	% 10eV

% classVsRel(radius, energy)
% waveFuncComparison(radius, energy)
% decompositionClassVsRel(radius, energy)

% fluxAndCurr(radius)
% relFluxAndCurr(radius)

%% Davidovic rods - star disk
function davidovicStarDisk
    % Create two-plot tile with a disk (radius, r=2m) rotating with
    % velocity, v=0 and v=1.5*10^8m/s.
    radius = 2;
    speed = [0 1.5*10^8];

    DavR = DavidovicRods(radius, speed);
    starDisk(DavR);
end

%% Paraboloidal deformation of the disk
function diskDef
    % Plot deformation of the disk surface under Lorenz contraction of the
    % circumference.
    
    % Visualli flat paraboloid
    semiAxes1 = [1000 1000];
    radius1 = 8;
    paraboloidHandle1 = Paraboloid(semiAxes1, 'centrePoint', [0 0 0],...
                                              'orientation', '-',...
                                              'meshDens', 1,...
                                              'rLim', radius1);
    paraboloid1 = getParaboloid(paraboloidHandle1).coordinates;
    circleHandle1 = Circle(radius1, 'centrePoint', [0 0 0]);
    circle1 = getCircle(circleHandle1).coordinates;
    radiusCoord1 = {[0 0] [-radius1 0] [0 0]}; % Radius coordinates to plot
    
    % Convex paraboloid
    semiAxes2 = [10 10];
    radius2 = 6;
    paraboloidHandle2 = Paraboloid(semiAxes2, 'centrePoint', [0 0 0],...
                                              'orientation', '-',...
                                              'meshDens', 1,...
                                              'rLim', radius2);
    paraboloid2 = getParaboloid(paraboloidHandle2).coordinates;
    circleHandle2 = Circle(radius2, 'centrePoint', [0 0 0]);
    circle2 = getCircle(circleHandle2).coordinates;
    % Radius coordinates to plot
    ry = -radius2:0.1:0;
    rx = zeros(1, length(ry));
    rz = -(ry.^2)./10 + (radius2^2)/10;
    radiusCoord2 = {rx ry rz};
    
    % Shift paraboloid to align its egdes with the circle - found
    % radius-shift dependency to be (radius2^2)/10.
    paraboloid2{3} = paraboloid2{3} + (radius2^2)/10;
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    % Define two tiles
    layout = defineTile(layout, 'axesNames', {'x [m]' 'y [m]' ''},...
                                'size', [-10 10 -10 10 -10 10],...
                                'legend', 'top',...
                                'title', "Static disc");
                            
    layout = defineTile(layout, 'axesNames', {'x [m]' 'y [m]' ''},...
                                'size', [-10 10 -10 10 -10 10],...
                                'legend', 'top',...
                                'title', "Disc rotating with relativistic velocity");
    
	% Add circles and paraboloids to the tiles
    layout = addPlot(layout, 1, circle1,        'color', 'b',...
                                                'lineWidth', 2,...
                                                'name', "Circumference");
                                            
	layout = addPlot(layout, 1, radiusCoord1,   'color', 'r',...
                                                'lineWidth', 1.5,...
                                                'name', "Radius");
                                            
    layout = addSurf(layout, 1, paraboloid1,    'faceColor', 'none',...
                                                'name', "Surface");
    
    layout = addPlot(layout, 2, circle2,        'color', 'b',...
                                                'lineWidth', 2,...
                                                'name', "Circumference");
                                            
    layout = addPlot(layout, 2, radiusCoord2,   'color', 'r',...
                                                'lineWidth', 1.5,...
                                                'name', "Radius");
                                            
    layout = addSurf(layout, 2, paraboloid2,    'faceColor', 'none',...
                                                'name', "Surface");
    
    drawLayout(layout);
end

%% Simple wavefunction plots
function wavefunction
    len = 1;          % Length of the box
    x = 0:0.01:len;   % X-domain
    % Wavefunction for sin and cos
    psiSinFunc = @(n) sqrt(2/len).*sin(n*pi*x/len);
    psiCosFunc = @(n) sqrt(2/len).*cos(n*pi*x/len);
    % Wavefunction domain cell arrays
    psiSin = {psiSinFunc(0) psiSinFunc(1) psiSinFunc(2) psiSinFunc(3) psiSinFunc(4)};
    psiCos = {psiCosFunc(0) psiCosFunc(1) psiCosFunc(2) psiCosFunc(3) psiCosFunc(4)};
    
    lim = findLimits(psiSinFunc(4)).*1.05;
    
    layout = Plot;
    layout = createLayout(layout, 5, 2);
    for i = 1:10
        layout = defineTile(layout, 'size', [0 len lim],...
                                    'legend', 'none',...
                                    'axesNames', {'' ''});
    end
    j = 1;
    for i = [1 3 5 7 9]
        layout = addPlot(layout, i, {x psiSin{j}}, 'lineWidth', 1,...
                                                   'color', 'k');
        j = j + 1;
    end
    j = 1;
    for i = [2 4 6 8 10]
        layout = addPlot(layout, i, {x, psiCos{j}}, 'lineWidth', 1,...
                                                    'color', 'k');
        j = j + 1;
    end
    drawLayout(layout)
end

%% Comparing different orders of perturbation approximation of the energy
function correctionOrders
    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    % 2st, 3nd and 4rd order approximations - correspond to the 2nd, 3rd
    % and 4th term of the power series
    energyApprox = {0 0 0};
    energyApprox{1} = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                  'order', 2).coordinates;
	energyApprox{2} = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                  'order', 3).coordinates;
    energyApprox{3} = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                  'order', 4).coordinates;
	name = ["1st order" "2nd order" "3rd order"];
    
    % Create layout
	layout = Plot;
    layout = createLayout(layout, 1, 2);
    layout = defineTile(layout, 'title', "Speed range" + newline + "from 0.01c to the 0.99c" + newline,...
                                'axesNames', {'Speed (v/c)' 'Energy [J]'},...
                                'legend', 'bottom-left');
    layout = defineTile(layout, 'title', "Speed range" + newline + "from 0.90c to the 0.99c" + newline,...
                                'axesNames', {'Speed (v/c)' 'Energy [J]'},...
                                'size', [0.90 0.99 -8000 7500],...
                                'legend', 'bottom-left');
    % Add plots to each tile
	for i = 1:3
        layout = addPlot(layout, 1, energyApprox{i}, 'name', name(i));
        layout = addPlot(layout, 2, energyApprox{i}, 'name', name(i));
    end
    
    drawLayout(layout);
    
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "2nd, 3rd and 3rd order" + newline + "approximation of the energy");
end

function correctionOrdersLog
    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    % 2st, 3nd and 4rd order approximations - correspond to the 2nd, 3rd
    % and 4th term of the power series
    exactValue = getEnergyApproximation(energyApproxHandle, 'exact').coordinates;
    energyApprox2 = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                'order', 2).coordinates;
	energyApprox3 = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                'order', 3).coordinates;
    energyApprox4 = getEnergyApproximation(energyApproxHandle, 'approximation', 'model', 'relativistic',...
                                                                                'order', 4).coordinates;
    
	% Plot data in a logarithmic scale applied on the y-axis
    semilogy(energyApprox2{1}, energyApprox2{2},...
             energyApprox3{1}, energyApprox3{2},...
             energyApprox4{1}, energyApprox4{2},...
             exactValue{1}, exactValue{2}, 'k');
	
	% Add properties
	lgd = legend("2nd order", "3rd order", "4th order", "Exact");
    lgd.Location = 'northwest';
    title("2nd, 3rd and 3rd order" + newline + "approximation of the energy." + newline + "Logarithmic scale applied on the y-axis");
    xlabel("Speed (v/c)");
    ylabel("Energy [J]");
    grid on;
end

function correctionDeviationLog
    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    % 2st, 3nd and 4rd order approximations - correspond to the 2nd, 3rd
    % and 4th terms deviations from the exact value of the energy
    exactValue = getEnergyApproximation(energyApproxHandle, 'exact').coordinates;
    energyDev2 = getEnergyApproximation(energyApproxHandle, 'deviation', 'model', 'relativistic',...
                                                                         'order', 2).coordinates;
	energyDev3 = getEnergyApproximation(energyApproxHandle, 'deviation', 'model', 'relativistic',...
                                                                         'order', 3).coordinates;
    energyDev4 = getEnergyApproximation(energyApproxHandle, 'deviation', 'model', 'relativistic',...
                                                                         'order', 4).coordinates;
                                                                       
		% Plot data in a logarithmic scale applied on the y-axis
    semilogy(energyDev2{1}, energyDev2{2},...
             energyDev3{1}, energyDev3{2},...
             energyDev4{1}, energyDev4{2},...
             exactValue{1}, exactValue{2}, 'k');
	
	% Add properties
	lgd = legend("2nd order", "3rd order", "4th order", "Exact");
    lgd.Location = 'northwest';
    title("2nd, 3rd and 3rd order approximation" + newline + "deviation from the energy value." + newline + "Logarithmic scale applied on the y-axis");
    xlabel("Speed (v/c)");
    ylabel("Energy deviation [J]");
    grid on;
end

function correctionDeviationRel2ndOrder
     energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    % 2st approximations deviation from the energy relative to its value
    energyDev2 = getEnergyApproximation(energyApproxHandle, 'deviationRel', 'model', 'relativistic',...
                                                                            'order', 2).coordinates;
	layout = Plot;
    layout = createLayout(layout, 1, 1);
    layout = defineTile(layout, 'title', "2nd order approximation deviation" + newline + "from the exact energy relative" + newline + "to the value of energy",...
                                'size', [0.4 0.99 0 45],...
                                'legend', 'none',...
                                'axesNames', {'Speed (v/c)' 'Relative energy deviation'});
    layout = addPlot(layout, 1, energyDev2);
    drawLayout(layout);
end

function classVsRelApprox
    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    relApprox = getEnergyApproximation(energyApproxHandle, 'deviationRel', 'model', 'relativistic',...
                                                                           'order', 2).coordinates;
    classApprox = getEnergyApproximation(energyApproxHandle, 'deviationRel', 'model', 'classical',...
                                                                             'order', 2).coordinates;
    
    layout = Plot;
    layout = createLayout(layout, 1, 1);
    layout = defineTile(layout, 'size', [0 1 0 3],...
                                'title', "2nd-order power series approximation" + newline +...
                                         "applied to relativistic and classical" + newline +...
                                         "models - relative deviation of the" + newline +...
                                         "approximation from the exact value of energy",...
                                'legend', 'top-left',...
                                'axesNames', {'Speed (v/c)' 'Relative energy deviation'});
    
    layout = addPlot(layout, 1, relApprox, 'name', "Relativistic");
    layout = addPlot(layout, 1, classApprox, 'name', "Classical",...
                                             'color', 'm');
    
    drawLayout(layout);
end

%% Wavefunction
function classVsRel(radius, energy)
    qFactor = 2000;
    classWave = {{} {} {}};
    relWave = {{} {} {}};
    E = [0.1*energy energy 10*energy];
    
    m = 9.1094*10.^(-31);                   % Electron rest mass
    c = 3*10^8;                             % Speed of light in vacuum
    Er = m*c^2;                             % Rest mass energy
    speed = @(E) c*sqrt(1-(Er/(E + Er))^2); % Speed of an electron
    
    for j = 1:3
        % Finding speed corresponding to the energy input
        speed(E(j))
        
        % Finding quantum numbers for given radius-energy input for classical
        % and relativistic model
        quantumN = QuantumN(radius, 'energy', E(j),...
                                    'relCorrection', false);
        classList = getTheList(quantumN)

        quantumN = QuantumN(radius, 'energy', E(j),...
                                    'relCorrection', true);
        relList = getTheList(quantumN)
    
        % Generating wavefunctions for given quantum numbers and superimposing
        % them
        class = {zeros(1, qFactor+1) zeros(1, qFactor+1)};
        rel = class;
    
        for i = 1:length(classList)
            wavefunctionHandle = Wavefunction(radius, classList(i), 'q', qFactor);
            coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;

            class{1} = class{1} + (1/sqrt((length(classList)))).*coordinates{1};
            class{2} = class{2} + (1/sqrt((length(classList)))).*coordinates{2};
        end
        for i = 1:length(relList)
            wavefunctionHandle = Wavefunction(radius, relList(i), 'q', qFactor);
            coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;

            rel{1} = rel{1} + (1/sqrt((length(relList)))).*coordinates{1};
            rel{2} = rel{2} + (1/sqrt((length(relList)))).*coordinates{2};
        end
        
        classWave{j} = class;
        relWave{j} = rel;
    end
    
    E = E./(1.602*10^-19);
    
    % Creating layout
    layout = Plot;
    layout = createLayout(layout, 3, 2);
    
    layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                'title', "Classical model." + newline + "E=" + string(E(1)) + "eV",...
                                'legend', 'none',...
                                'size', [-2200 4800 -3200 3200]);
    layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                'title', "Relativistic model." + newline + "E=" + string(E(1)) + "eV",...
                                'legend', 'none',...
                                'size', [-2200 4800 -3200 3200]);
	for j = 1:2
        layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                    'title', "E=" + string(E(j+1)) + "eV",...
                                    'legend', 'none',...
                                    'size', [-2200 4800 -3200 3200]);
        layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                    'title', "E=" + string(E(j+1)) + "eV",...
                                    'legend', 'none',...
                                    'size', [-2200 4800 -3200 3200]);
    end

    for j = 1:3
        layout = addPlot(layout, (2*j-1), classWave{j}, 'color', '#666666');
        layout = addPlot(layout, (2*j), relWave{j}, 'color', '#666666');
    end
	drawLayout(layout);
	
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Wavefunction plotted in the x-y plane." + newline +...
             "Amplitude as an extension of the radial direction.");
end

function waveFuncComparison(radius, energy)
    % Finding quantum numbers for given radius-energy input
    quantumN = QuantumN(radius, 'energy', energy,...
                                'relCorrection', true);
    list = getTheList(quantumN);
    
    qFactor = 2000; % Number of steps for the function to take when
                    % generating coordinate arrays
                    
	circleHandle = Circle(radius, 'centrePoint', [0 0 0]);
    circle = getCircle(circleHandle).coordinates;
    
                    
    % Preallocating arrays for wavefunction components
    wavefunction = {zeros(1, qFactor+1) zeros(1, qFactor+1)};
    wave = {zeros(1, qFactor+1) zeros(1, qFactor+1) zeros(1, qFactor+1)};
    
    for i = 1:length(list)
        wavefunctionHandle = Wavefunction(radius, list(i), 'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;
        
        wavefunction{1} = wavefunction{1} + (1/sqrt((length(list)))).*coordinates{1};
        wavefunction{2} = wavefunction{2} + (1/sqrt((length(list)))).*coordinates{2};
    end
    
    for i = 1:length(list)
        waveHandle = Wave(radius, list(i), 'q', qFactor);
        coordinates = getWave(waveHandle, 'arithmeticType', 'cos').coordinates;
        
        wave{3} = wave{3} + (1/sqrt((length(list)))).*coordinates{3};
    end
    
    wave{1} = coordinates{1};
    wave{2} = coordinates{2};
    
    layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                'title', "Wavefunction plotted in the x-y plane." + newline +...
                                         "Amplitude as an extension of the radial" + newline +...
                                         "direction.");
    layout = defineTile(layout, 'axesNames', {'x [m]' 'y [m]' 'Amplitude [m^-^1^/^2]'},...
                                'title', "Wavefunction plotted in the three-dimensional" + newline +...
                                         "space. Electron path in the x-y plane and" + newline +...
                                         "amplitude in the z direction.");
    
    layout = addPlot(layout, 1, wavefunction, 'color', '#666666',...
                                              'name', "Wavefunction");
    layout = addPlot(layout, 1, {circle{1} circle{2}}, 'lineSpec', '-.r',...
                                                       'lineWidth', 2,...
                                                       'name', "Electron's path");
    
    layout = addPlot(layout, 2, wave, 'color', '#666666',...
                                      'name', "Wavefunction");
    layout = addPlot(layout, 2, circle, 'lineSpec', '-.r',...
                                        'lineWidth', 2,...
                                        'name', "Electron's path");
    
    drawLayout(layout);
    
end

function decompositionClassVsRel(radius, energy)
    qFactor = 2000;
    classWave = {zeros(1, qFactor+1) zeros(1, qFactor+1) zeros(1, qFactor+1)};
    relWave = classWave;
    energy = 10*energy;

    % Finding quantum numbers for given radius-energy input for classical
    % and relativistic model
    quantumN = QuantumN(radius, 'energy', energy,...
                                'relCorrection', false);
    classList = getTheList(quantumN)

    quantumN = QuantumN(radius, 'energy', energy,...
                                'relCorrection', true);
    relList = getTheList(quantumN)
    
    % Generating wavefunctions for given quantum numbers and superimposing
    % them
    for i = 1:length(classList)
        wavefunctionHandle = Wavefunction(radius, classList(i), 'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;

        classWave{1} = classWave{1} + (1/sqrt((length(classList)))).*coordinates{1};
        classWave{2} = classWave{2} + (1/sqrt((length(classList)))).*coordinates{2};
        
        waveHandle = Wave(radius, classList(i), 'q', qFactor);
        coordinates = getWave(waveHandle, 'arithmeticType', 'cos').coordinates;
        
        classWave{3} = classWave{3} + (1/sqrt((length(classList)))).*coordinates{3};
    end
    for i = 1:length(relList)
        wavefunctionHandle = Wavefunction(radius, relList(i), 'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;

        relWave{1} = relWave{1} + (1/sqrt((length(relList)))).*coordinates{1};
        relWave{2} = relWave{2} + (1/sqrt((length(relList)))).*coordinates{2};
        
        waveHandle = Wave(radius, relList(i), 'q', qFactor);
        coordinates = getWave(waveHandle, 'arithmeticType', 'cos').coordinates;
        
        relWave{3} = relWave{3} + (1/sqrt((length(relList)))).*coordinates{3};
    end
    
    % Creating layout
    layout = Plot;
    layout = createLayout(layout, 2, 2);
    
    layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                'title', "Classical model.",...
                                'size', [-1100 5300 -3600 3600],...
                                'legend', 'none');
    layout = defineTile(layout, 'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m] + Amplitude [m^-^1^/^2]'},...
                                'title', "Relativistic model.",...
                                'size', [-1100 5300 -3600 3600],...
                                'legend', 'none');
    layout = defineTile(layout, 'axesNames', {'Angle [rad]' 'Amplitude [m^-^1^/^2]'},...
                                'title', 'none',...
                                'size', [-inf inf -3900 5300],...
                                'legend', 'none');
    layout = defineTile(layout, 'axesNames', {'Angle [rad]' 'Amplitude [m^-^1^/^2]'},...
                                'title', 'none',...
                                'size', [-inf inf -3900 5300],...
                                'legend', 'none');
                            
    layout = addPlot(layout, 1, {classWave{1} classWave{2}}, 'color', '#666666');
    layout = addPlot(layout, 2, {relWave{1} relWave{2}}, 'color', '#666666');
    layout = addPlot(layout, 3, {0:2*pi/qFactor:2*pi classWave{3}}, 'color', [0.4660, 0.6740, 0.1880]);
    layout = addPlot(layout, 4, {0:2*pi/qFactor:2*pi relWave{3}}, 'color', [0.4660, 0.6740, 0.1880]);
    
    drawLayout(layout);
                            
	% Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Wavefunction in x-y plane (top) and its amplitude" + newline +...
             "vs angle (bottom). Classical (left) vs relativistic" + newline +...
             "model (right). High energy input, E=100eV.");
end

%% Current density & magnetic flux
function fluxAndCurr(radius)
    currDensHandle = CurrentDensity(radius, 'relCorrection', true);
    relCurrDens = getCurrentDensity(currDensHandle).coordinates;
    maxVel = getCurrentDensity(currDensHandle).maxVel;
    
    currDensHandle = CurrentDensity(radius, 'relCorrection', false);
    classCurrDens = getCurrentDensity(currDensHandle).coordinates;
    
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', true);
    relFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', false);
    classFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    % Change units keV, kT, mA/rad
    relFluxDens = {relFluxDens{1}./10^3 relFluxDens{2}./10^3};
    classFluxDens = {classFluxDens{1}./10^3 classFluxDens{2}./10^3};
    relCurrDens = {relCurrDens{1}./10^3 relCurrDens{2}.*10^3};
    classCurrDens = {classCurrDens{1}./10^3 classCurrDens{2}.*10^3};
    
    
    % Creating layout
    layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    layout = defineTile(layout, 'axesNames', {'Energy [keV]' 'Current density [mA\cdotrad^-^1]'},...
                                'title', "Current density",...
                                'legend', 'top-left');
    layout = defineTile(layout, 'axesNames', {'Energy [keV]' 'Flux density [kT]'},...
                                'title', "Magnetic flux density",...
                                'legend', 'top-left');
                            
    layout = addPlot(layout, 1, relCurrDens, 'lineSpec', 'o',...
                                             'name', "Applied correction to energy",...
                                             'color', [0.6350, 0.0780, 0.1840]);
    layout = addPlot(layout, 1, classCurrDens, 'lineSpec', 'o',...
                                               'name', "Classical model",...
                                               'color', [0.3010, 0.7450, 0.9330]);
	layout = addPlot(layout, 2, relFluxDens, 'lineSpec', 'o',...
                                             'name', "Applied correction to energy",...
                                             'color', [0.4940, 0.1840, 0.5560]);
    layout = addPlot(layout, 2, classFluxDens, 'lineSpec', 'o',...
                                               'name', "Classical model",...
                                               'color', [0.4660, 0.6740, 0.1880]); 
	
    drawLayout(layout);
    
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Current density and magnetic flux density with and" + newline + "without the perturbation correction applied.");
end

function relFluxAndCurr(radius)
    currDensHandle = CurrentDensity(radius, 'relCorrection', true);
    relCurrDens = getCurrentDensity(currDensHandle).coordinates;
    
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', true);
    relFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    % Change units keV, kT, uA/rad
    relFluxDens = {relFluxDens{1}./10^3 relFluxDens{2}./(150*10^3)};
    relCurrDens = {relCurrDens{1}./10^3 relCurrDens{2}.*10^3};
    
    % Creating layout
    layout = Plot;
    layout = createLayout(layout, 1, 1);
    
    layout = defineTile(layout, 'axesNames', {"Energy [keV]" "Current density [mA\cdotrad^-^1]" + newline + "Flux density [150\cdotkT]"},...
                                'title', "Current density and magnetic flux density with and" + newline + "without the perturbation correction applied.",...
                                'legend', 'top-left');
                            
    layout = addPlot(layout, 1, relCurrDens, 'lineSpec', 'o',...
                                             'name', "Current density",...
                                             'color', [0.6350, 0.0780, 0.1840]);
	layout = addPlot(layout, 1, relFluxDens, 'lineSpec', 'o',...
                                             'name', "Magnetic flux density",...
                                             'color', [0.4940, 0.1840, 0.5560]);
	
    drawLayout(layout);
end
