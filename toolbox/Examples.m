format compact

% circleEx1
% circleEx2
% spiralEx1
% spiralEx2
% paraboloidEx1
% paraboloidEx2
% wavefunctionEx
% plotEx
% quantumNEx1
% quantumNEx2
% waveEx
% gifEx
% davidovicRodsEx
% energyApproximationEx
% currentDensityEx
magneticFluxEx

%% CIRCLE (ex. 1)
function circleEx1
    % Create a circle with radius 5 and default CIRCLE input
    
    radius = 5;
    
    circleHandle = Circle(radius);
    circle = getCircle(circleHandle);
    
    plot(circle.coordinates{1}, circle.coordinates{2});
    axis(circle.size.*2);
    
    view(2)
    xlabel('x')
    ylabel('y')
    grid on
end

%% CIRCLE (ex. 2)
function circleEx2
    % Create a circle in 3D (fixed x coordinate) with radius 3.14, centred
    % at [0, 3.14, 3.14] - it is plotted on the y-z surface at x = 0.
    
    radius = 3.14;
    
    circleHandle = Circle(radius, 'centrePoint', [0 3.14 3.14],...
                                  'fixedCoordinate', 'x');
    circle = getCircle(circleHandle);
    
    plot3(circle.coordinates{1}, circle.coordinates{2}, circle.coordinates{3});
    axis(circle.size);
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end

%% SPIRAL (ex. 1)
function spiralEx1
    % Create a spiral with radii range from 3 to 5 an default SPIRAL input
    
    rMin = 3;
    rMax = 5;
    
    spiralHandle = Spiral(rMin, rMax);
	spiral = getSpiral(spiralHandle);
    
	plot(spiral.coordinates{1}, spiral.coordinates{2});
	axis(spiral.size);
    
    view(2)
    xlabel('x')
    ylabel('y')
    grid on
end

%% SPIRAL (ex. 2)
function spiralEx2
    % Create a spiral in 3D (fixed x coordinate) with radius ranging from
    % 6.28 to 12.56 centred at [0, 3.14, 3.14] - it is plotted on the y-z 
    % surface at x = 0 - and 4 laps.
    
    rMin = 6.28;
    rMax = 12.56;
    
    spiralHandle = Spiral(rMin, rMax, 'centrePoint', [0 3.14 3.14],...
                                      'fixedCoordinate', 'x',...
                                      'noOfLaps', 4,...
                                      'q', 180);
	spiral = getSpiral(spiralHandle);
    
	plot3(spiral.coordinates{1}, spiral.coordinates{2}, spiral.coordinates{3});
	axis(spiral.size);
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end

%% PARABOLOID (ex. 1)
function paraboloidEx1
    % Create a paraboloid centred at x = 1, y = 3 with default class
    % parameters and mesh density 0.5.

    semiAxes = [1 3];
    
    paraboloidHandle = Paraboloid(semiAxes, 'meshDens', 0.5);
    paraboloid = getParaboloid(paraboloidHandle);
    
    surf(paraboloid.coordinates{1}, paraboloid.coordinates{2}, paraboloid.coordinates{3});
	axis(paraboloid.size)
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end

%% PARABOLOID (ex. 2)
function paraboloidEx2
    % Create reversed hyperbolic paraboloid centred at x = 1, y = 1.

    semiAxes = [1 1];
    
    paraboloidHandle = Paraboloid(semiAxes, 'centrePoint', [3 3 1],...
                                            'orientation', '-',...
                                            'type', 'hyperbolic');
    paraboloid = getParaboloid(paraboloidHandle);
    
    surf(paraboloid.coordinates{1}, paraboloid.coordinates{2}, paraboloid.coordinates{3});
	axis(paraboloid.size)
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end

%% WAVEFUNCTION (ex.)
function wavefunctionEx
    % Plot 'sin' component of a wavefunction for the time, t=0 and quantum
    % state described by n=6 for an electron travelling around the circle
    % with radius, r=0.1

    radius = 0.1;
    quantumN = 6;
    time = 1;
    
    wavefunctionHandle = Wavefunction(radius, quantumN);
	wavefunction = getWavefunc(wavefunctionHandle, time);
    
	plot3(wavefunction.coordinates{1}, wavefunction.coordinates{2}, wavefunction.coordinates{3});
	axis(wavefunction.size)
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end

%% PLOT (ex.)
function plotEx
    % Create two circles
    radius1 = 5;
    radius2 = 3;
    
    circleHandle1 = Circle(radius1, 'centrePoint', [0 0]);
    circle1 = getCircle(circleHandle1);
    
    circleHandle2 = Circle(radius2, 'centrePoint', [1 1 1]);
    circle2 = getCircle(circleHandle2);
    
    % Create a spiral
    rMin = 1;
    rMax = 5;
    
	spiralHandle = Spiral(rMin, rMax, 'noOfLaps', 2);
	spiral = getSpiral(spiralHandle);
    
    % Create a paraboloid
    semiAxes = [2 2];
    
	paraboloidHandle = Paraboloid(semiAxes, 'meshDens', 2);
	paraboloid = getParaboloid(paraboloidHandle);
    
    % Define m-by-n layout (two tiles)
    m = 1;
    n = 2;
    
    layout = Plot;
	layout = createLayout(layout, m, n);
    
    % Define tiles
    layout = defineTile(layout, 'title',        "Circle + Spiral",...
                                'axesNames',    {'x', 'y'},...
                                'size',         circle1.size.*2);
	layout = defineTile(layout, 'title',        "Circle + Paraboloid",...
                                'axesNames',    {'x', 'y', 'z'},...
                                'size',         circle2.size.*2);
                            
    % Add plots/surfaces to tiles
	layout = addPlot(layout, 1, circle1.coordinates, 'lineSpec',    "--r",...
                                                     'name',        "Circle");
	layout = addPlot(layout, 2, circle2.coordinates, 'color',       "#555555",...
                                                     'name',        "Circle");
	layout = addPlot(layout, 1, spiral.coordinates, 'color',        'b',...
                                                    'name',         "Spiral");
	layout = addSurf(layout, 2, paraboloid.coordinates, 'edgeColor',    [0.8 0.5 0.1],...
                                                        'faceColor',    'none',...
                                                        'faceAlpha',    0,...
                                                        'name',         "Paraboloid");
    
    % Draw the layout
    drawLayout(layout);
end

%% QUANTUMN (ex. 1)
function quantumNEx1
    % An electron moves with speed 15000m/s around a circular path of
    % radius 3000nm 
    radius = 0.000003;

    quantumN = QuantumN(radius, 'speed', 15000);
	list = getTheList(quantumN);
    
    disp(list);
end

%% QUANTUMN (ex. 2)
function quantumNEx2
    % Determine the list of allowed quantum numbers when an electron with
    % energy of 3*10^(-19) joules. The radius is 0.1mm, relativisitic
    % correction is applied.
    radius = 0.0001;

    quantumN = QuantumN(radius, 'energy', 3*10^(-19),...
                                'relCorrection', true);
    list = getTheList(quantumN);

    % Generating wavefunction for each quantum number
    time = 0;
    
    % Superimposing wavefunction in each x, y, z direction
    qFactor = 2000; % Number of steps for the function to take when
                    % generating coordinate arrays
    
    sumx0 = zeros(1, qFactor+1);
    sumy0 = zeros(1, qFactor+1);
    sumz0 = zeros(1, qFactor+1);
    
    for i = 1:length(list)
        wavefunctionHandle = Wavefunction(radius, list(i), 'q', qFactor);
        coordinates = getWavefunc(wavefunctionHandle, time).coordinates;
        limits = getWavefunc(wavefunctionHandle, time).size.*2;
        
        sumx0 = sumx0 + (1/sqrt((length(list)))).*coordinates{1};
        sumy0 = sumy0 + (1/sqrt((length(list)))).*coordinates{2};
        sumz0 = sumz0 + (1/sqrt((length(list)))).*coordinates{3};
    end

    % Define m-by-n layout (two tiles)
    m = 1;
    n = 2;

    layout = Plot;
    layout = createLayout(layout, m, n);

    % Define tiles
    layout = defineTile(layout, 'title', "Superimposed wavefunction",...
                                'axesNames', {'x', 'y', 'Imaginary axis'},...
                                'size', limits,...
                                'legend', 'none');
    layout = defineTile(layout, 'title', "Decomposed wavefunction",...
                                'axesNames', {'Angle (deg)', 'Spatial coordinate'});

    % Add plots to tiles
    pltArray1 = {sumx0 sumy0 sumz0};
    pltArray2x = {1:qFactor+1 sumx0};
    pltArray2y = {1:qFactor+1 sumy0};
    pltArray2z = {1:qFactor+1 sumz0};
    
    layout = addPlot(layout, 1, pltArray1, 'color', 'k',...
                                           'name', "Wavefunction");
    layout = addPlot(layout, 2, pltArray2x, 'color', 'c',...
                                            'name', "X-component");
    layout = addPlot(layout, 2, pltArray2y, 'color', 'm',...
                                            'name', "Y-component");
    layout = addPlot(layout, 2, pltArray2z, 'color', 'g',...
                                            'name', "Imaginary component");

    % Draw the layout
    drawLayout(layout);
end

%% WAVE (ex.)
function waveEx
    % Plot a wavefunction of an electron which energy is described by
    % n=8' It is travelling around the circular path of radius, r=0.1m.
    % Desired type of the wavefunction is cosine at time, t=0
    
    radius = 0.1;
    quantumN = 8;
    time = 0;
    
    % Create a circle - wave horizontal axis
    circleHandle = Circle(radius, 'centrePoint', [0 0 0]);
    circle = getCircle(circleHandle);

    % Create a wave
    waveHandle = Wave(radius, quantumN);
    wave = getWave(waveHandle, time, 'arithmeticType', 'cos');
    
    plot3(wave.coordinates{1}, wave.coordinates{2}, wave.coordinates{3}, 'k', circle.coordinates{1}, circle.coordinates{2}, circle.coordinates{3}, '--');
    axis(wave.size)
    
    view(3)
    title("Superimposed wavefunction")
    legend("Wavefunction", "Trajectory")
    xlabel('x')
    ylabel('y')
    zlabel('Amplitude')
    grid on
end

%% GIF (ex.)
function gifEx
    
end

%% DAVIDOVICRODS (ex.)
function davidovicRodsEx
    radius = 2;
    speed = 0*10^8:0.5*10^8:2.5*10^8;

    DavR = DavidovicRods(radius, speed);
    rodLength(DavR);
end

%% ENERGYAPPROXIMATION (ex.)
function energyApproximationEx
    % Compare 2nd and 3rd order power series (classical) approximation to
    % the electrons energy moving with the speed higher than 0.1c and
    % smaller than 0.9c. The speed resolution is 0.01.
    vMin = 0.1;
    vMax = 0.9;
    vRes = 0.01;

    enApproxHandle = EnergyApproximation(vMin, vMax, vRes);
    enApprox1 = getEnergyApproximation(enApproxHandle, 'approximation', 'model', 'classical',...
                                                                        'order', 2);
    enApprox2 = getEnergyApproximation(enApproxHandle, 'approximation', 'model', 'classical',...
                                                                        'order', 3);
    
    plot(enApprox1.coordinates{1}, enApprox1.coordinates{2},...
         enApprox2.coordinates{1}, enApprox2.coordinates{2});
    
    view(2)
    xlabel('Relative speed')
    ylabel('Energy')
    grid on
end

%% CURRENTDENSITY (ex.)
function currentDensityEx
    % Plot change of the current density across the allowed speed range
    radius = 0.00000003;
    
    currentDensityHandle = CurrentDensity(radius, 'relCorrection', true);
	currentDensity = getCurrentDensity(currentDensityHandle);
    
    plot(currentDensity.coordinates{1}, currentDensity.coordinates{2}, 'o');
end

%% MAGNETICFLUX (ex.)
function magneticFluxEx
    % Plot change in the quantised magnetic flux density across the allowed
    % range of speed values
    radius = 0.00000003;
    
    magneticFluxHandle = MagneticFlux(radius, 'relCorrection', true);
    magneticFlux = getMagneticFlux(magneticFluxHandle);
    
    plot(magneticFlux.coordinates{1}, magneticFlux.coordinates{2}, 'o');
end
