format compact

% davidovicStarDisk
% diskDef
vincentWavefunction

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
function vincentWavefunction
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
