%% RESULTS
% The following code was written to for the "Investigation of the Ehrenfest
% Paradox for a Particle on a Ring Model" 3rd year project. Later edited to
% fit the page width.
%
% Patryk Jesionka, 2020
%%

format compact

% davidovicStarDisk

% diskDef

% wavefunction

% correctionOrders
% correctionDeviationRel2ndOrder
% classVsRelApprox

radius = 35*10^(-9);        % 35nm
energy = 10*1.602*10^-19;	% 10eV

% classVsRel(radius, energy)
% waveFuncComparison(radius, energy)
% decompositionClassVsRel(radius, energy)

% fluxAndCurr(radius)
% relFluxAndCurr(radius)
% qualityFactor(100)
% qualityFactor(1000)
% qualityFactor(10000)

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
    layout = defineTile(layout,...
        'axesNames', {'x [m]' 'y [m]' ''},...
        'size', [-10 10 -10 10 -10 10],...
        'legend', 'top',...
        'title', "Static disc");
                            
    layout = defineTile(layout,...
        'axesNames', {'x [m]' 'y [m]' ''},...
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
    % Plot sine and cosine wavefunctions described in the Vincent's
    % publication in the 5-by-2 tiled layout.

    len = 1;          % Length of the box
    x = 0:0.01:len;   % X-domain
    
    % Wavefunction for sin and cos
    psiSinFunc = @(n) sqrt(2/len).*sin(n*pi*x/len);
    psiCosFunc = @(n) sqrt(2/len).*cos(n*pi*x/len);
    
    % Wavefunction domain cell arrays
    psiSin = {psiSinFunc(0) psiSinFunc(1)...
        psiSinFunc(2) psiSinFunc(3) psiSinFunc(4)};
    psiCos = {psiCosFunc(0) psiCosFunc(1)...
        psiCosFunc(2) psiCosFunc(3) psiCosFunc(4)};
    
    % Find appropriate axes sizes
    lim = findLimits(psiSinFunc(4)).*1.05;
    
    % Create a layout
    layout = Plot;
    layout = createLayout(layout, 5, 2);
    
    % Define tiles
    for i = 1:10
        layout = defineTile(layout,...
            'size', [0 len lim],...
            'legend', 'none',...
            'axesNames', {'' ''});
    end
    
    % Add sine component wavefunctions to the tile
    j = 1;
    for i = [1 3 5 7 9]
        layout = addPlot(layout, i, {x psiSin{j}},...
            'lineWidth', 1,...
            'color', 'k');
        j = j + 1;
    end
    
    % Add cosine component wavefunctions to the tile
    j = 1;
    for i = [2 4 6 8 10]
        layout = addPlot(layout, i, {x, psiCos{j}},...
            'lineWidth', 1,...
            'color', 'k');
        j = j + 1;
    end
    
    % Draw layout
    drawLayout(layout)
end

%% Comparing different orders of perturbation approximation of the energy
function correctionOrders
    % Define handle for the EnergyApproximation object
    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    
    % Exact value of the energy
    exactValue = getEnergyApproximation(energyApproxHandle,...
        'exact').coordinates;
    
    % 2st, 3nd and 4rd order approximations
    energyApprox = {0 0 0};
    energyApprox{1} = getEnergyApproximation(energyApproxHandle,...
        'approximation',...
        'model', 'relativistic',...
        'order', 2).coordinates;
	energyApprox{2} = getEnergyApproximation(energyApproxHandle,...
        'approximation',...
        'model', 'relativistic',...
        'order', 3).coordinates;
    energyApprox{3} = getEnergyApproximation(energyApproxHandle,...
        'approximation',...
        'model', 'relativistic',...
        'order', 4).coordinates;
	
    % Define names of each element
    name = ["2nd order" "3rd order" "4th order" "Exact"];
    
    % Create layout
	layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'title', "Speed range" + newline + "from 0.01c to the 0.99c"...
        + newline,...
        'axesNames', {'Speed (v/c)' 'Energy [J]'},...
        'size', [0.01 0.99 0 1.2*10^-13],...
        'legend', 'top-left');
    layout = defineTile(layout,...
        'title', "Speed range" + newline + "from 0.01c to the 0.78c"...
        + newline,...
        'axesNames', {'Speed (v/c)' 'Energy [J]'},...
        'size', [0.4 0.7 0 0.4*10^-13],...
        'legend', 'top-left');
    
    % Add plots to each tile
	for i = 1:3
        layout = addPlot(layout, 1, energyApprox{i}, 'name', name(i));
        layout = addPlot(layout, 2, energyApprox{i}, 'name', name(i));
    end
    layout = addPlot(layout, 1, exactValue, 'name', name(4),...
                                            'color', 'k');
    layout = addPlot(layout, 2, exactValue, 'name', name(4),...
                                            'color', 'k');
    
    drawLayout(layout);
    
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "2nd, 3rd and 3rd order" + newline + "approximation of the"...
        + " energy vs" + newline + "exact value of the energy");
end

function correctionDeviationRel2ndOrder
    % Plot 2nd order approximation deviation from the exact value of the
    % energy relative to it

    c = 2.998*10^8; % Speed of light in the vacuum

    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    
    % 2st approximations deviation from the energy relative to its value
    energyDev2 = getEnergyApproximation(energyApproxHandle,...
        'deviationRel',...
        'model', 'relativistic',...
        'order', 2).coordinates;
    
    % Create layout
	layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'size', [0.01 0.99 0 20],...
        'legend', 'none',...
        'axesNames', {'Speed (v/c)' 'Relative energy deviation'});
    layout = defineTile(layout,...
        'size', [0.7 0.99 0 1.1],...
        'legend', 'none',...
        'axesNames', {'Speed (v/c)' 'Relative energy deviation'});
    
    % Add plots
    layout = addPlot(layout, 1, energyDev2);
    layout = addPlot(layout, 2, energyDev2);
    layout = addPlot(layout, 2, {[0.894 0.894] [0 1]},...
        'color', '#666666',...
        'lineStyle', '-.');
	layout = addPlot(layout, 2, {[0.928 0.928] [0 0.1729]},...
        'color', '#666666',...
        'lineStyle', '-.');
    
    drawLayout(layout);
    
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "2nd order approximation deviation" + newline +...
             "from the exact energy relative" + newline +...
             "to the value of energy");
         
    % Add description of the local minima and maxima
	text(0.894, 1.03, changeNotationType(0.89*c, 's') + "ms^-^1" +...
        newline, 'HorizontalAlignment', 'center');
    text(0.928, 0.1729, " " + changeNotationType(0.93*c, 's') +...
        "ms^-^1", 'HorizontalAlignment', 'left');
end

function classVsRelApprox
    % Compare relative deviation to its value with classical momentum input

    energyApproxHandle = EnergyApproximation(0.01, 0.99, 'step', 0.002);
    relApprox = getEnergyApproximation(energyApproxHandle,...
        'deviationRel',...
        'model', 'relativistic',...
        'order', 2).coordinates;
    classApprox = getEnergyApproximation(energyApproxHandle,...
        'deviationRel',...
        'model', 'classical',...
        'order', 2).coordinates;
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 1, 1);
    
    % Define tiles
    layout = defineTile(layout,...
        'size', [0 1 0 3],...
        'title', "2nd-order power series approximation" + newline +...
            "applied to relativistic and classical" + newline +...
            "models - relative deviation of the" + newline +...
            "approximation from the exact value of energy",...
        'legend', 'top-left',...
        'axesNames', {'Speed (v/c)' 'Relative energy deviation'});
    
	% Add plots
    layout = addPlot(layout, 1, relApprox, 'name', "Relativistic");
    layout = addPlot(layout, 1, classApprox, 'name', "Classical",...
                                             'color', 'm');
    
    drawLayout(layout);
end

%% Wavefunction
function classVsRel(radius, energy)
    % Accuracy of the plot
    qFactor = 2000;
    
    % Prelocated cell arrays
    classWave = {{} {} {}};
    relWave = {{} {} {}};
    
    % Energies to be fed into the functions
    E = [0.1*energy energy 10*energy];
    
    % Define constants
    m = 9.1094*10.^(-31);                   % Electron rest mass
    c = 3*10^8;                             % Speed of light in vacuum
    Er = m*c^2;                             % Rest mass energy
    
    speed = @(E) c*sqrt(1-(Er/(E + Er))^2); % Speed of an electron
    
    % Loop through each energy value from the matrix
    for j = 1:3
        % Finding speed corresponding to the energy input
        speed(E(j))
        
        % Finding quantum numbers for given radius-energy input for
        % classical and relativistic model
        quantumN = QuantumN(radius, 'energy', E(j),...
                                    'relCorrection', false);
        classList = getTheList(quantumN)

        quantumN = QuantumN(radius, 'energy', E(j),...
                                    'relCorrection', true);
        relList = getTheList(quantumN)
    
        % Generating wavefunctions for given quantum numbers and
        % superimposing them
        class = {zeros(1, qFactor+1) zeros(1, qFactor+1)};
        rel = class;
    
        % Loop through each electrionic state for classical and
        % relativistic model
        for i = 1:length(classList)
            wavefunctionHandle = Wavefunction(radius, classList(i),...
                'q', qFactor);
            coordinates = getWavefunc(wavefunctionHandle,...
                'arithmeticType', 'cos').coordinates;

            class{1} = class{1} + (1/sqrt((length(classList)))).*coordinates{1};
            class{2} = class{2} + (1/sqrt((length(classList)))).*coordinates{2};
        end
        for i = 1:length(relList)
            wavefunctionHandle = Wavefunction(radius, relList(i), 'q', qFactor);
            coordinates = getWavefunc(wavefunctionHandle, 'arithmeticType', 'cos').coordinates;

            rel{1} = rel{1} + (1/sqrt((length(relList)))).*coordinates{1};
            rel{2} = rel{2} + (1/sqrt((length(relList)))).*coordinates{2};
        end
        
        % Add values to cells
        classWave{j} = class;
        relWave{j} = rel;
    end
    
    % Units, J -> eV
    E = E./(1.602*10^-19);
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 3, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
        ' + Amplitude [m^-^1^/^2]'},...
        'title', "Classical model." + newline + "E=" +...
        string(E(1)) + "eV",...
        'legend', 'none',...
        'size', [-2200 4800 -3200 3200]);
    layout = defineTile(layout,...
        'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
        ' + Amplitude [m^-^1^/^2]'},...
        'title', "Relativistic model." + newline + "E=" +...
        string(E(1)) + "eV",...
        'legend', 'none',...
        'size', [-2200 4800 -3200 3200]);
                            
	for j = 1:2
        layout = defineTile(layout,...
            'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
            ' + Amplitude [m^-^1^/^2]'},...
            'title', "E=" + string(E(j+1)) + "eV",...
            'legend', 'none',...
            'size', [-2200 4800 -3200 3200]);
        layout = defineTile(layout,...
            'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
            ' + Amplitude [m^-^1^/^2]'},...
            'title', "E=" + string(E(j+1)) + "eV",...
            'legend', 'none',...
            'size', [-2200 4800 -3200 3200]);
    end

    % Add plots
    for j = 1:3
        layout = addPlot(layout, (2*j-1), classWave{j},...
            'color', '#666666');
        layout = addPlot(layout, (2*j), relWave{j},...
            'color', '#666666');
    end
    
	drawLayout(layout);
	
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Wavefunction plotted in the x-y plane." + newline +...
             "Amplitude as an extension of the radial direction.");
end

function waveFuncComparison(radius, energy)
    % Compare two different ways of plotting wavefunction

    % Finding quantum numbers for given radius-energy input
    quantumN = QuantumN(radius, 'energy', energy,...
                                'relCorrection', true);
    list = getTheList(quantumN);
    
    qFactor = 2000;
    
    % Electron's trajectory
	circleHandle = Circle(radius, 'centrePoint', [0 0 0]);
    circle = getCircle(circleHandle).coordinates;
    
    % Preallocating arrays for wavefunction components
    wavefunction = {zeros(1, qFactor+1) zeros(1, qFactor+1)};
    wave = {zeros(1, qFactor+1) zeros(1, qFactor+1) zeros(1, qFactor+1)};
    
    % Get and superimpose wavefunctions
    for i = 1:length(list)
        wavefunctionHandle = Wavefunction(radius, list(i), 'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle,...
            'arithmeticType', 'cos').coordinates;
        
        wavefunction{1} = wavefunction{1} + (1/sqrt((length(list))))...
            .*coordinates{1};
        wavefunction{2} = wavefunction{2} + (1/sqrt((length(list))))...
            .*coordinates{2};
    end
    
    for i = 1:length(list)
        waveHandle = Wave(radius, list(i), 'q', qFactor);
        coordinates = getWave(waveHandle,...
            'arithmeticType', 'cos').coordinates;
        
        wave{3} = wave{3} + (1/sqrt((length(list)))).*coordinates{3};
    end
    
    % x and y coordinates for the WAVE remain unchanged
    wave{1} = coordinates{1};
    wave{2} = coordinates{2};
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
        ' + Amplitude [m^-^1^/^2]'},...
        'title', "Wavefunction plotted in the x-y plane." + newline +...
        "Amplitude as an extension of the radial" + newline +...
        "direction.");
    layout = defineTile(layout,...
        'axesNames', {'x [m]' 'y [m]' 'Amplitude [m^-^1^/^2]'},...
        'title', "Wavefunction plotted in the three-dimensional" +...
        newline + "space. Electron path in the x-y plane" + newline +...
        "and amplitude in the z direction.");
    
    % Add plots
    layout = addPlot(layout, 1, wavefunction,...
        'color', '#666666',...
        'name', "Wavefunction");
    layout = addPlot(layout, 1, {circle{1} circle{2}},...
        'lineSpec', '-.r',...
        'lineWidth', 2,...
        'name', "Electron's path");
    
    layout = addPlot(layout, 2, wave,...
        'color', '#666666',...
        'name', "Wavefunction");
    layout = addPlot(layout, 2, circle,...
        'lineSpec', '-.r',...
        'lineWidth', 2,...
        'name', "Electron's path");
    
    drawLayout(layout);
end

function decompositionClassVsRel(radius, energy)
    % Show the wavefunction fro high-energy input and its amplitude plotted
    % against the angle

    qFactor = 2000;
    classWave = {zeros(1, qFactor+1)...
                 zeros(1, qFactor+1)...
                 zeros(1, qFactor+1)};
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
        wavefunctionHandle = Wavefunction(radius, classList(i),...
            'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle,...
            'arithmeticType', 'cos').coordinates;
        
        % Superposition
        classWave{1} = classWave{1} + (1/sqrt((length(classList))))...
            .*coordinates{1};
        classWave{2} = classWave{2} + (1/sqrt((length(classList))))...
            .*coordinates{2};
        
        waveHandle = Wave(radius, classList(i), 'q', qFactor);
        coordinates = getWave(waveHandle,...
            'arithmeticType', 'cos').coordinates;
        
        % Superposition
        classWave{3} = classWave{3} + (1/sqrt((length(classList))))...
            .*coordinates{3};
    end
    for i = 1:length(relList)
        wavefunctionHandle = Wavefunction(radius, relList(i),...
            'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle,...
            'arithmeticType', 'cos').coordinates;

        % Superposition
        relWave{1} = relWave{1} + (1/sqrt((length(relList))))...
            .*coordinates{1};
        relWave{2} = relWave{2} + (1/sqrt((length(relList))))...
            .*coordinates{2};
        
        waveHandle = Wave(radius, relList(i), 'q', qFactor);
        coordinates = getWave(waveHandle,...
            'arithmeticType', 'cos').coordinates;
        
        % Superposition
        relWave{3} = relWave{3} + (1/sqrt((length(relList))))...
            .*coordinates{3};
    end
    
    % Creat layout
    layout = Plot;
    layout = createLayout(layout, 2, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
        ' + Amplitude [m^-^1^/^2]'},...
        'title', "Classical model.",...
        'size', [-1100 5300 -3600 3600],...
        'legend', 'none');
    layout = defineTile(layout,...
        'axesNames', {'x [m] + Amplitude [m^-^1^/^2]' 'y [m]' +...
        ' + Amplitude [m^-^1^/^2]'},...
        'title', "Relativistic model.",...
        'size', [-1100 5300 -3600 3600],...
        'legend', 'none');
    layout = defineTile(layout,...
        'axesNames', {'Angle [rad]' 'Amplitude [m^-^1^/^2]'},...
        'title', 'none',...
        'size', [-inf inf -3900 5300],...
        'legend', 'none');
    layout = defineTile(layout,...
        'axesNames', {'Angle [rad]' 'Amplitude [m^-^1^/^2]'},...
        'title', 'none',...
        'size', [-inf inf -3900 5300],...
        'legend', 'none');
	
    % Add plots
    layout = addPlot(layout, 1, {classWave{1} classWave{2}},...
        'color', '#666666');
    layout = addPlot(layout, 2, {relWave{1} relWave{2}},...
    'color', '#666666');
    layout = addPlot(layout, 3, {0:2*pi/qFactor:2*pi classWave{3}},...
        'color', [0.4660, 0.6740, 0.1880]);
    layout = addPlot(layout, 4, {0:2*pi/qFactor:2*pi relWave{3}},...
        'color', [0.4660, 0.6740, 0.1880]);
    
    drawLayout(layout);
                            
	% Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Wavefunction in x-y plane (top) and its amplitude" +...
             newline + "vs angle (bottom). Classical (left) vs" +...
             " relativistic" + newline + "model (right). High energy" +...
             " input, E=100eV.");
end

%% Current density & magnetic flux
function fluxAndCurr(radius)
    % Plot current density and magnetic flux density - classical and
    % relativistic cases
    
    % Relativistic current density
    currDensHandle = CurrentDensity(radius, 'relCorrection', true);
    relCurrDens = getCurrentDensity(currDensHandle).coordinates;
    
    % Output data
    maxVel = getCurrentDensity(currDensHandle).maxVel
    maxRelCurr = findLimits(relCurrDens{2})
    maxEnergy = findLimits(relCurrDens{1})
    
    % Classical current density
    currDensHandle = CurrentDensity(radius, 'relCorrection', false);
    classCurrDens = getCurrentDensity(currDensHandle).coordinates;
    
    % Output data
    maxClassCurr = findLimits(classCurrDens{2})
    
    % Relativistic magnetic flux density
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', true);
    relFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    % Output data
    maxRelFlux = findLimits(relFluxDens{2})
    
    % Classical magnetic flux density
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', false);
    classFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    % Output data
    maxClassFlux = findLimits(classFluxDens{2})
    
    % Change units keV, kT, mA/rad
    relFluxDens = {relFluxDens{1}./10^3 relFluxDens{2}./10^3};
    classFluxDens = {classFluxDens{1}./10^3 classFluxDens{2}./10^3};
    relCurrDens = {relCurrDens{1}./10^3 relCurrDens{2}.*10^3};
    classCurrDens = {classCurrDens{1}./10^3 classCurrDens{2}.*10^3};
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 1, 2);
    
    % Define tiles
    layout = defineTile(layout,...
        'axesNames', {'Energy [keV]'...
        'Current density [mA\cdotrad^-^1]'},...
        'title', "Current density",...
        'legend', 'top-left');
    layout = defineTile(layout,...
        'axesNames', {'Energy [keV]' 'Flux density [kT]'},...
        'title', "Magnetic flux density",...
        'legend', 'top-left');
                           
	% Add plots
    layout = addPlot(layout, 1, relCurrDens,...
        'lineSpec', 'o',...
        'name', "Applied correction to energy",...
        'color', [0.6350, 0.0780, 0.1840]);
    layout = addPlot(layout, 1, classCurrDens,...
        'lineSpec', 'o',...
        'name', "Classical model",...
        'color', [0.3010, 0.7450, 0.9330]);
	layout = addPlot(layout, 2, relFluxDens,...
        'lineSpec', 'o',...
        'name', "Applied correction to energy",...
        'color', [0.4940, 0.1840, 0.5560]);
    layout = addPlot(layout, 2, classFluxDens,...
        'lineSpec', 'o',...
        'name', "Classical model",...
        'color', [0.4660, 0.6740, 0.1880]); 
	
    drawLayout(layout);
    
    % Add title to the layout
    t = findobj('type', 'tiledlayout');
    title(t, "Current density and magnetic flux density with and" + newline + "without the perturbation correction applied.");
end

function relFluxAndCurr(radius)
    % Plot relativisitc current and magnetic flux densities

    % Current density
    currDensHandle = CurrentDensity(radius, 'relCorrection', true);
    relCurrDens = getCurrentDensity(currDensHandle).coordinates;
    
    % Magnetic flux
    fluxDensHandle = MagneticFlux(radius, 'relCorrection', true);
    relFluxDens = getMagneticFlux(fluxDensHandle).coordinates;
    
    % Change units keV, kT, uA/rad
    relFluxDens = {relFluxDens{1}./10^3 relFluxDens{2}./(150*10^3)};
    relCurrDens = {relCurrDens{1}./10^3 relCurrDens{2}.*10^3};
    
    % Create layout
    layout = Plot;
    layout = createLayout(layout, 1, 1);
    
    % Define tile
    layout = defineTile(layout,...
        'axesNames', {"Energy [keV]" "Current density [mA\cdotrad^-^1]"...
        + newline + "Flux density [150\cdotkT]"},...
        'title', "Current density and magnetic flux density with"...
        + newline + "the perturbation correction applied.",...
        'legend', 'top-left');
                            
	% Add plots
    layout = addPlot(layout, 1, relCurrDens,...
        'lineSpec', 'o',...
        'name', "Current density",...
        'color', [0.6350, 0.0780, 0.1840]);
	layout = addPlot(layout, 1, relFluxDens,...
        'lineSpec', 'o',...
        'name', "Magnetic flux density",...
        'color', [0.4940, 0.1840, 0.5560]);
	
    drawLayout(layout);
end

%% Quality factor
function qualityFactor(q)
    % Determine the list of allowed quantum numbers when an electron with
    % energy of 3*10^(-19) joules. The radius is 0.1mm, relativisitic
    % correction is applied. Then plot the case for given quality factor.
    
    radius = 0.0001;

    quantumN = QuantumN(radius, 'energy', 3*10^(-19),...
                                'relCorrection', true);
    list = getTheList(quantumN);
    
    % Superimposing wavefunction in each x, y, z direction
    sumx0 = zeros(1, q+1);
    sumy0 = zeros(1, q+1);
    
    for i = 1:length(list)
        wavefunctionHandle = Wavefunction(radius, list(i), 'q', q);
        coordinates = getWavefunc(wavefunctionHandle,...
            'arithmeticType', 'cos').coordinates;
        
        sumx0 = sumx0 + (1/sqrt((length(list)))).*coordinates{1};
        sumy0 = sumy0 + (1/sqrt((length(list)))).*coordinates{2};
    end

    % Define 1-by-1 layout (one tile)
    layout = Plot;
    layout = createLayout(layout, 1, 2);

    % Define tiles
    layout = defineTile(layout,...
        'title', "q=" + string(q),...
        'axesNames', {'Steps', 'Amplitude [m^-^1^/^2]'});
	layout = defineTile(layout,...
        'title', "Enlarged, q=" + string(q),...
        'axesNames', {'Steps', 'Amplitude [m^-^1^/^2]'},...
        'size', [(4/10)*q (5/10)*q -30 30]);

    % Add plots to the tile
    pltArray2x = {1:q+1 sumx0};
    pltArray2y = {1:q+1 sumy0};
    for i = 1:2
        layout = addPlot(layout, i, pltArray2x, 'color', 'g',...
                                                'name', "X-component");
        layout = addPlot(layout, i, pltArray2y, 'color', 'm',...
                                                'name', "Y-component");
    end

    % Draw the layout
    drawLayout(layout);
end
