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

%% CIRCLE (ex. 1)
function circleEx1
    % Create a circle in 2D (x-y plane) with radius 5, centred at [0, 0]
    
    radius = 5;
    centrePoint = [0 0];
    fixedCoordinate = 'z';
    
    circleHandle = Circle(radius, centrePoint, fixedCoordinate);
    circle = getCircle(circleHandle);
    
    plot(circle.coordinates{1}, circle.coordinates{2});
    axis(circle.size);
    
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
    centrePoint = [0 3.14 3.14];
    fixedCoordinate = 'x';
    
    circleHandle = Circle(radius, centrePoint, fixedCoordinate);
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
    % Create a spiral in 2D (x-y plane) with radii range from 3 to 5, 
    % centred at [0, 0] and 1 lap.
    
    rMin = 3;
    rMax = 5;
    centrePoint = [0 0];
    fixedCoordinate = 'z';
    noOfLaps = 1;
    
    spiralHandle = Spiral(rMin, rMax, centrePoint, fixedCoordinate, noOfLaps);
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
    centrePoint = [0 3.14 3.14];
    fixedCoordinate = 'x';
    noOfLaps = 4;
    
    spiralHandle = Spiral(rMin, rMax, centrePoint, fixedCoordinate, noOfLaps);
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
    centrePoint = [0 0 0];
    semiAxes = [1 3];
    orientation = '+';
    meshDens = 5;
    
    paraboloidHandle = Paraboloid(centrePoint, semiAxes, orientation, meshDens);
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
    centrePoint = [3 3 1];
    semiAxes = [1 -1];
    orientation = '+';
    meshDens = 25;
    
    paraboloidHandle = Paraboloid(centrePoint, semiAxes, orientation, meshDens);
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
    radius = 0.1;
    quantumN = 6;
    time = 0;
    arithmeticType = 'sin';
    amplitudeAxes = 'xy';
    
    
    wavefunctionHandle = Wavefunction(radius, quantumN);
	wavefunction = getWavefunc(wavefunctionHandle, time, arithmeticType, amplitudeAxes);
    
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
    centrePoint1 = [0 0];
    fixedCoordinate1 = 'z';
    radius2 = 3;
    centrePoint2 = [1 1 1];
    fixedCoordinate2 = 'z';
    
    circleHandle1 = Circle(radius1, centrePoint1, fixedCoordinate1);
    circle1 = getCircle(circleHandle1);
    circleHandle2 = Circle(radius2, centrePoint2, fixedCoordinate2);
    circle2 = getCircle(circleHandle2);
    
    % Create a spiral
    rMin = 1;
    rMax = 5;
    centrePoint = [0 0];
    fixedCoordinate = 'z';
    noOfLaps = 2;
    
	spiralHandle = Spiral(rMin, rMax, centrePoint, fixedCoordinate, noOfLaps);
	spiral = getSpiral(spiralHandle);
    
    % Create a paraboloid
    centrePoint = [0 0 0];
    semiAxes = [2 2];
    orientation = '+';
    meshDens = 2;
    
	paraboloidHandle = Paraboloid(centrePoint, semiAxes, orientation, meshDens);
	paraboloid = getParaboloid(paraboloidHandle);
    
    % Define m-by-n layout (two tiles)
    m = 1;
    n = 2;
    
    layout = Plot;
	layout = createLayout(layout, m, n);
    
    % Define tiles
    title1 = "Circle + Spiral";
    axesNames1 = {'x', 'y'};
    size1 = circle1.size;
    title2 = "Circle + Paraboloid";
    axesNames2 = {'x', 'y', 'z'};
    size2 = circle2.size;
    
	layout = defineTile(layout, title1, axesNames1, size1);
	layout = defineTile(layout, title2, axesNames2, size2);
    
    % Add plots/surfaces to tiles
	layout = addPlot(layout, 1, circle1.coordinates, '--', 'r', 0.5, "Circle");
	layout = addPlot(layout, 2, circle2.coordinates, '-', 'r', 0.5, "Circle");
	layout = addPlot(layout, 1, spiral.coordinates, '-', 'b', 0.5, "Spiral");
	layout = addSurf(layout, 2, paraboloid.coordinates, 'k', '-', 'none', 0, "Paraboloid");
    
    % Draw the layout
    drawLayout(layout);
end

%% QUANTUMN (ex. 1)
function quantumNEx1
    radius = 0.000003;
    linearSpeed = 15000;

    quantumN = QuantumN(radius, linearSpeed);
	list = getTheList(quantumN);
    
    disp(list);
end

%% QUANTUMN (ex. 2)
function quantumNEx2
    % Determine the list of allowed quantum numbers
    linearSpeed = 15000;
    radius = 0.000003;

    quantumN = QuantumN(linearSpeed, radius);
    list = getTheList(quantumN);

    % Generating wavefunction for each quantum number
    time = 0;
    arythmeticType = 'sin';
    wavefunction = [];

    for i = size(list)
        wavefunctionHandle = Wavefunction(radius, list(1));
        wavefunction = [wavefunction, getWavefunc(wavefunctionHandle, time, arythmeticType).coordinates];
    end

    % Superimposing wavefunction in each x, y, z direction
    sumx = zeros(1, 201);
    sumy = zeros(1, 201);
    sumz = zeros(1, 201);
    for i = 1:3:size(list)
        sumx = sumx + wavefunction{i};
        sumy = sumy + wavefunction{i+1};
        sumz = sumz + wavefunction{i+2};
    end

    % Define m-by-n layout (two tiles)
    m = 1;
    n = 2;

    layout = Plot;
    layout = createLayout(layout, m, n);

    % Define tiles
    title1 = "Superimposed wavefunction";
    axesNames1 = {'x', 'y', 'Imaginary axis'};
    size1 = [0 0 0 0 0 0]; % forces: axis auto
    title2 = "Decomposed wavefunction";
    axesNames2 = {' ', 'Spatial coordinate'};
    size2 = [0 201 -400 400];

    layout = defineTile(layout, title1, axesNames1, size1);
    layout = defineTile(layout, title2, axesNames2, size2);

    % Add plots to tiles
    layout = addPlot(layout, 1, {sumx sumy sumz}, '-', 'k', 0.5, "Wavefunction");
    layout = addPlot(layout, 2, {1:201 sumx}, '-', 'c', 0.5, "X-component");
    layout = addPlot(layout, 2, {1:201 sumy}, '-', 'm', 0.5, "Y-component");
    layout = addPlot(layout, 2, {1:201 sumz}, '-', 'g', 0.5, "Imaginary component");

    % Draw the layout
    drawLayout(layout);
end

%% WAVE (ex.)
function waveEx
    radius = 0.1;
    quantumN = 8;
    time = 0;
    arithmeticType = 'cos';
    
    % Create a circle - wave horizontal axis
    circleHandle = Circle(radius, [0 0 0], 'z');
    circle = getCircle(circleHandle);

    % Create a wave
    waveHandle = Wave(radius, quantumN);
    wave = getWave(waveHandle, time, arithmeticType);
    
    plot3(wave.coordinates{1}, wave.coordinates{2}, wave.coordinates{3}, circle.coordinates{1}, circle.coordinates{2}, circle.coordinates{3}, '--');
    axis(wave.size)
    
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
end