%% PLOT
% Class making an extensive use of other toolbox classes and plot, plot3,
% surf and tiledLayout Matlab functionalities. It generates defined-size
% layout containing tiles which may display 2D or 3D plots or surfaces
% with Matlab styling properties. 
%
%   obj = Plot
%       Constructor.
%
%           Output:
%       'obj': object of the class.
%
%   obj = createLayout(obj, m, n)
%       Uses built-in (tiledLayout() function) validation and creates m by
%       n tiled layout, saves number of tiles and assignes False flags
%       (since tiles are not defined yet) to each of them.
%       
%           Input:
%       'obj':          object of the class.
%       'm', 'n':       (required), positive integer values, size of a
%                       tiledLayout.
%           Output:
%       'obj':          object of the class.
%
%   obj = defineTile(obj)
%   obj = defineTile(obj, Name, Value)
%       Validates input, assigns data to the tile structure, then adds it
%       to the tiles array and updates definedTiles flag.
%       
%           Input:
%       'obj':          object of the class.
%       'title':        'none' (default), string of characters, title of a
%                       tile.
%       'axesNames':    {'x' 'y'} (default), array of 2 or 3 char or string
%                       type entries, names of each of the axis of the
%                       tile.
%       'size':         [-inf inf -inf inf] (default), array of 4 or 6
%                       numerical entries, size of the tile.
%       'legend':       'top-right' (default), location of the legend. Used
%                       more intuitive notation, e.g. 'top' instead of
%                       'north'.
%           Output:
%       'obj':          object of the class.
%
%   obj = addPlot(obj, tileNo, pltArray)
%   obj = addPlot(obj, tileNo, pltArray, Name, Value)
%       Validates input, adds plot to the tile existing in the tiles array.
%
%           Input:
%       'obj':          object of the class.
%       'tileNo':       (required), integer value between 1 and n*m,
%                       dictating to which tile the data passed to the
%                       function must be assigned.
%       'pltArray':     (required), 2 or 3 dimensional nested array with
%                       coordinates of the plot.
%       'lineSpec':     '-' (default), LineSpec property of plot/plot3.
%       'lineStyle':    'default' (default), LineStyle property for
%                       plot/plot3.
%       'color':        'auto' (default), Color property of plot/plot3.
%       'lineWidth':    0.5 (default), LineWidth property of plot/plot3.
%       'marker':       'none' (default), Marker property of plot/plot3.
%       'markerColor':  'auto' (default), MarkerColor property of
%                       plot/plot3.
%       'name':         'none' (default), name of the plot to be displayed
%                       in the legend.
%           Output:
%       'obj':          object of the class.
%
%   obj = addSurf(obj, tileNo, srfArray)
%   obj = addSurf(obj, tileNo, srfArray, Name, Value)
%       Validates input, adds plot to the tile existing in the tiles array.
%
%           Input:
%       'obj':          object of the class.
%       'tileNo':       (required), integer value between 1 and n*m,
%                       dictating to which tile the data passed to the
%                       function must be assigned.
%       'srfArray':     (required), 3 dimensional array containing mesh
%                       (matrices) of the surface.
%       'edgeColor': 	'auto' (default), EdgeColor property of surf.
%       'lineStyle':    '-' (default), LineStyle property of surf.
%       'lineWidth':    0.5 (default), LineWidth property of surf.
%       'faceColor':    'auto' (default), FaceColor property of surf.
%       'faceAlpha':    1 (default), FaceAlpha property of surf.
%       'name':         'none' (default), name of the surface to be
%                       displayed in the legend.
%           Output:
%       'obj':          object of the class.
%
%   drawLayout(obj)
%       Loops through each tile, uses drawTile function.
%
%           Input:
%       'obj':          object of the class.
%
%   handles = drawTile(obj, tile)
%       (private method) loops through plots/surfaces assigned to the tile,
%       uses drawPlot and drawSurf functions.
%
%           Input:
%       'obj':          object of the class.
%       'tile':         (required), tile structure.
%           Output:
%       'handles':      handles of plot/surf objects plotted on the tile.
%
%   handle = drawPlot(obj, object, dimensions)
%       (private method) uses plot() and plot3() Matlab functions and adds
%       style properties to the graph.
%
%           Input:
%       'obj':          object of the class.
%       'object':       (required), plot structure containing coordinates
%                       array and style properties.
%       'dimensions':   (required), numerical value (2 or 3), dictates
%                       dimensions of the plot.
%           Output:
%       'handle':       handle to the plot drawn.
%
%   handle = drawSurf(obj, object)
%       (private method) uses surf() Matlab function and adds style
%       properties to the graph.
%
%           Input:
%       'obj':          object of the class.
%       'object':       (required), surface structure containing mesh array
%                       and style properties.
%           Output:
%       'handle':       handle to the surface drawn.
%
%   flag = validColorFunc(~, color)
%       (private method) validate if the string or char array are proper 
%       color expressions.
%
%           Input:
%       'obj':          object of the class.
%       'color':        (required), variable with color expression to
%                       validate.
%           Output:
%       'flag':         true if the color has valid form, false otherwise.
%
%	flag = isDefault(~, property, expression)
%       Check if the property is set to the default value.
%
%           Input:
%       'property':     (required), the variable to verify.
%       'expression':   (required), an expression the property's default
%                       value is defined by, e.g. if property is 'none' by
%                       default, then expression = 'none'.
%           Output:
%       'flag':         true if the property = expression, false otherwise.
%
%   Limitations:
%       The class can only be fed with certain type of structure -
%       struct('coordinates', [], 'size', []). Therefore, it is recommended
%       to use other classes of the EhrenfestToolbox to generate data.
%       The class cannot generate animations (gifs).
%
%   Examples:
%       Function plotEx from Examples.m
%
%   Updates:
%       01/04/2020: Addeded support for Marker line specification.
%       02/05/2020: Added input parser and support for additional plot/surf
%           properties.
%       03/05/2020: Added support for excluding plot/surf names from a
%           legend.
%
%   Use:
%       It is recommended to use other classes of the EhrenfestToolbox to 
%       generate data, though, it can be fed with any arrays of appropriate
%       dimensionality - e.g. if the array with coordinates is
%       two-dimensional then the array with axes size must have 4 entries.
%           Call in order:
%       1. createLayout
%       2. defineTile (as many times as needed)
%       3. addPlot, addSurf (as many times as needed)
%       4. drawLayout 
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, GIF, DAVIDOVICRODS, CHANGENOTATIONTYPE,
%       FINDLIMITS, CURRENTDENSITY, MAGNETICFLUX
%
%   Patryk Jesionka, 2019
%%

classdef Plot
    properties
        layout              % Tiled layout handle
        
        noOfTiles = 1       % Number of tiles available - used to validate
                            % if user specified existing tile_no
                            
        definedTiles = {}   % Contains information whether the certain tile
                            % was defined or not
        
        tiles = {}          % Array containing structures for every tile
        
        % Structure containing plots data
        plt = struct('coordinatesArray', [],...
                     'lineSpec', [],...
                     'lineStyle', [],...
                     'color', [],...
                     'lineWidth', [],...
                     'marker', [],...
                     'markerColor', [],... 
                     'name', [],...
                     'noOfFrames', [],...
                     'type', 'plot')
                 
        srf = struct('surfaceArray', [],...
                     'edgeColor', [],...
                     'lineStyle', [],...
                     'lineWidth', [],...
                     'faceColor', [],...
                     'faceAlpha', [],...
                     'name', [],...
                     'noOfFrames', [],...
                     'type', 'surface')
    end
    methods (Access = private)
        function handle = drawPlot(obj, object, dimensions)
            % Creating plot object (first assigning coordinates for the
            % sake of explicity)
            x = object.coordinatesArray{1};
            y = object.coordinatesArray{2};
            if dimensions == 3
                z = object.coordinatesArray{3};
                p = plot3(x, y, z, object.lineSpec);
            else
                p = plot(x, y, object.lineSpec);
            end
            
            % Setting properties
            if ~isDefault(obj, object.lineStyle, 'default')
                p.LineStyle = object.lineStyle;
            end
            if ~isDefault(obj, object.color, 'auto')
                p.Color = object.color;
            end
            
            p.LineWidth = object.lineWidth;
            
            if ~isDefault(obj, object.marker, 'none')
                p.Marker = object.marker;
            end
            if ~isDefault(obj, object.markerColor, 'auto')
                p.MarkerEdgeColor = object.markerColor;
            end
            if ~isDefault(obj, object.name, 'none')
                p.DisplayName = object.name;
            end
            
            % Return plot handle
            handle = p;
        end
        
        function handle = drawSurf(obj, object)
            % Creating surface object (first assigning coordinates for the
            % sake of explicity)
            X = object.surfaceArray{1};
            Y = object.surfaceArray{2};
            Z = object.surfaceArray{3};
            s = surf(X, Y, Z);
            
            % Setting properties
            if ~isDefault(obj, object.edgeColor, 'auto')
                s.EdgeColor = object.edgeColor;
            end
            
            s.LineStyle = object.lineStyle;
            s.LineWidth = object.lineWidth;
            
            if ~isDefault(obj, object.faceColor, 'auto')
                s.FaceColor = object.faceColor;
            end
            
            s.FaceAlpha = object.faceAlpha;
            
            if ~isDefault(obj, object.name, 'none')
                s.DisplayName = object.name;
            end
            
            % Return surf handle
            handle = s;
        end
        
        function handles = drawTile(obj, tile)
            handles = [];
            
            % Looping through plots
            for j = 1:length(tile.plots)
                % If it is a 2d plot
                if tile.dimensions == 2
                    handles = [handles, drawPlot(obj, tile.plots{j}, 2)];
                elseif tile.dimensions == 3 % if it is three-dimensional
                    if tile.plots{j}.type == "plot"
                        handles = [handles, drawPlot(obj, tile.plots{j}, 3)];
                    elseif tile.plots{j}.type == "surface"
                        handles = [handles, drawSurf(obj, tile.plots{j})];
                    end
                end
            end
        end
        
        function flag = validColorFunc(~, color)
            % Function checks if the string or char array are proper for
            % MATLAB color expressions.
            
            % Define memberships sets
            colorShortNameSet = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'w'};
            colorNameSet = {'red' 'green' 'blue'...
                            'cyan' 'magenta' 'yellow'...
                            'black' 'white' 'none'};
            hexSet = ['0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'A'...
                      'a' 'B' 'b' 'C' 'c' 'D' 'd' 'E' 'e' 'F' 'f'];
            
            % Anonymous validation functions for colors
            validRGBTriplet = @(x) isrow(x) && eq(length(x), 3) && all(ge(x, 0)) && all(ge(1, x)) && all(isreal(x));
            validHexCode = @(x) ischar(x) && eq(length(x), 7) && eq(x(1), '#') && all(ismember(x(2:7), hexSet));
            validColor = @(x) ismember(x, colorNameSet) || ismember(x, colorShortNameSet) || validHexCode(x) || eq(x, 'auto');
            
            if isstring(color) || ischar(color)
                % Conversion from string to char
                if isstring(color)
                    color = char(color);
                end
                % Validation of the char
                if validColor(color)
                    flag = true;
                else
                    flag = false;
                end
            else
                % Validation if the value is RGB triplet
                if validRGBTriplet(color)
                    flag = true;
                else
                    flag = false;
                end
            end
        end
        
        function flag = isDefault(~, property, expression)
            % Check if the property has a value set by default. Some cases
            % concern values different than char array 'default',
            % therefore, the expression property specifying the value that
            % the property should be compared with.
            
            if isstring(property) || ischar(property)
                if ismember(property, {expression})
                    flag = true;
                else
                    flag = false;
                end
            else
                flag = false;
            end
        end 
    end
    methods (Access = public)
        function obj = createLayout(obj, m, n)
            % m - "height", n - "width"

            % Validation functions
            validMN = @(x) gt(x, 0) && isreal(x) && isnumeric(x) && isfinite(x) && isscalar(x) && eq(x, floor(x));
            
            if validMN(m) && validMN(n)
                % Creates m by n layout; then, specifies number of tiles
                obj.layout = tiledlayout(m, n);
                obj.noOfTiles = m*n;

                % Iterates through the number of tiles and adds false 
                % values since tiles are not defined yet
                for i = 1:obj.noOfTiles
                    obj.definedTiles{i} = false;
                end
            end
        end
        
        function obj = defineTile(obj, varargin)
            % Define default values
            defaultTitle = 'none';
            defaultAxesNames = {'x' 'y'};
            defaultSize = [-inf inf -inf inf];
            defaultLegend = 'top-right';
            
            % Define membership sets - legendLocationKeySet defines more
            % understandable options in terms of location of the legend,
            % thus, they are to be specified by the user, then translated
            % into the MATLAB values - stored in the legendLocationValueSet
            legendLocationKeySet = {'top' 'bottom' 'right' 'left'...
                                    'top-right' 'top-left'...
                                    'bottom-right' 'bottom-left'...
                                    'top-out' 'bottom-out' 'right-out'...
                                    'left-out' 'top-right-out'...
                                    'top-left-out' 'bottom-right-out'...
                                    'bottom-left-out' 'best' 'best-out'...
                                    'none'};
            legendLocationValueSet = {'north' 'south' 'east' 'west'...
                                      'northeast' 'northwest'...
                                      'southeast' 'southwest'...
                                      'northoutside' 'southoutside'...
                                      'eastoutside' 'westoutside'...
                                      'northeastoutside'...
                                      'northwestoutside'...
                                      'southeastoutside'...
                                      'southwestoutside' 'best'...
                                      'bestoutside' 'none'};
            
            % Validation functions
            validTitle = @(x) isstring(x) || ischar(x);
            validAxesNames = @(x) iscell(x) && isrow(x) && (eq(length(x), 2) || eq(length(x), 3));
            validSize = @(x) all(isnumeric(x)) && isrow(x) && (eq(length(x), 4) || eq(length(x), 6));
            validLegend = @(x) ismember(x, legendLocationKeySet);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;

            % Adding arguments
            addParameter(p, 'title', defaultTitle, validTitle);
            addParameter(p, 'axesNames', defaultAxesNames, validAxesNames);
            addParameter(p, 'size', defaultSize, validSize);
            addParameter(p, 'legend', defaultLegend, validLegend);

            parse(p, varargin{:});

            % Extract variables from the parser
            title = p.Results.title;
            axesNames = p.Results.axesNames;
            size = p.Results.size;
            legend = p.Results.legend;
            
            % Adjusting default length of the size vector if the plot is
            % meant to be three-dimensional instead of two.
            if length(size) == length(defaultSize)
                if all(eq(size, defaultSize))
                    % If size was not specified
                    if length(axesNames) == 3
                        % If the plot is meant to be three-dimensional
                        size = [-inf inf -inf inf -inf inf];
                    end
                end
            end   
            
            % Additional validation for the axesNames 
            if ~(2*length(axesNames) == length(size))
                error("axesNames must be half the length of size.");
            end
            for i = 1:length(axesNames)
                if ~(ischar(axesNames{i}) || isstring(axesNames{i}) || isnumeric(axesNames{i}))
                    error("axesName must be char, string or numeric type!")
                end
            end
            
            % Convert legend into a member of legendLocationValueSet
            m = containers.Map(legendLocationKeySet, legendLocationValueSet);
            m = values(m, { legend });
            legend = m{1};
            
            % Check if there is a space for another tile
            if ~(length(obj.tiles)+1 <= obj.noOfTiles)
                error("There is no space for another tile in the layout!");
            end
            
            % Assigning data to the tile structure
            tile.title = title;
            tile.axesNames = axesNames;
            tile.size = size;
            tile.legend = legend;
            tile.dimensions = length(axesNames);
            tile.plots = {};
            tile.handles = [];
            
            % Adds structure to the tile array
            obj.tiles{length(obj.tiles)+1} = tile;
            
            % Set the first false flag in definedTiles to true
            for i = 1:obj.noOfTiles
                if ~obj.definedTiles{i}
                    obj.definedTiles{i} = true;
                end
            end
        end
        
        function obj = addPlot(obj, tileNo, pltArray, varargin)
            % Define default values
            defaultLineSpec = '-';
            defaultLineStyle = 'default';
            defaultColor = 'auto';
            defaultLineWidth = 0.5;
            defaultMarker = 'none';
            defaultMarkerColor = 'auto';
            defaultName = 'none';
            
            % Define memberships sets
            lineStyleSet = {'-' '--' ':' '-.' ''};
            markerSet = {'o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' ''};
            colorShortNameSet = {'r' 'g' 'b' 'c' 'm' 'y' 'k' 'w' ''};
            
            % LineSpec membership set requires all possible combinations of
            % lineStyleSet, markerSet and colorShortNameSet elements
            n = length(lineStyleSet)*length(markerSet)*length(colorShortNameSet);
            lineSpecSet = cell(1, n);

            m = 1;
            for i = 1:length(lineStyleSet)
                for j = 1:length(markerSet)
                    for k = 1:length(colorShortNameSet)
                        lineSpecSet{m} = [lineStyleSet{i} markerSet{j} colorShortNameSet{k}];
                        m = m + 1;
                    end
                end
            end

            % Validation functions
            validTileNo = @(x) ge(x, 1) && isreal(x) && isnumeric(x) && eq(floor(x), x) && ge(obj.noOfTiles, x);
            validPltArray = @(x) iscell(x) && isrow(x) && ge(length(x), 2);
            validLineSpec = @(x) ismember(x, lineSpecSet);
            validLineStyle = @(x) ismember(x, [lineStyleSet 'default']);
            validColor = @(x) validColorFunc(obj, x);
            validLineWidth = @(x) ge(x, 0) && isreal(x) && isnumeric(x) && isscalar(x) && isfinite(x);
            validMarker = @(x) ismember(x, [markerSet 'none']);
            validMarkerColor = @(x) validColorFunc(obj, x);
            validName = @(x) isstring(x) || ischar(x);
            
            % Validation function for each inner cell of pltArray -
            % assuming that the number length(pltArray) is greater than
            % three (potential gif frames).
            validGifFrame = @(x1, x2) isrow(x1) && isrow(x2) && iscell(x1) && iscell(x2) && eq(length(x1), length(x2)) && (eq(length(x1), 2) || eq(length(x1), 3));
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;

            % Adding arguments
            addRequired(p, 'tileNo', validTileNo);
            addRequired(p, 'pltArray', validPltArray);
            addParameter(p, 'lineSpec', defaultLineSpec, validLineSpec);
            addParameter(p, 'lineStyle', defaultLineStyle, validLineStyle);
            addParameter(p, 'color', defaultColor, validColor);
            addParameter(p, 'lineWidth', defaultLineWidth, validLineWidth);
            addParameter(p, 'marker', defaultMarker, validMarker);
            addParameter(p, 'markerColor', defaultMarkerColor, validMarkerColor);
            addParameter(p, 'name', defaultName, validName);

            parse(p, tileNo, pltArray, varargin{:});
            
            % Additional validation tileNo and pltArray length
            if ~obj.definedTiles{p.Results.tileNo}
                error("This tile is not defined!");
            elseif ~(obj.tiles{p.Results.tileNo}.dimensions == length(p.Results.pltArray))
                error("Tile has different dimensions than the plot you want to add!");
            end
            
            % Additional validation for the pltArray - check if the input
            % was specified for the gif animation
            if gt(length(p.Results.pltArray), 3)
                % It is potential gif input
                for i = 2:length(p.Results.pltArray)
                    % Check for validity of each cell according to
                    % validGifFrame function
                    if ~validGifFrame(p.Results.pltArray{i-1}, p.Results.pltArray{i})
                        error("Input pltArray is neither proper plot nor gif input.");
                    end
                end
                % Check if the number of frames is consistent with
                % other plots
                if gt(length(obj.tiles{tileNo}), 0) && ~eq(length(p.Results.pltArray), obj.tiles{tileNo}.plots{1}.noOfFrames)
                    error("Number of frames is too small/large for this tile.");
                else
                    % Set number of frames
                    obj.plt.noOfFrames = length(p.Results.pltArray);
                end
            else
                % Input is for a plot
                obj.plt.noOfFrames = 1;
            end

            % Extract variables from the parser
            tileNo = p.Results.tileNo;
            obj.plt.coordinatesArray = p.Results.pltArray;
            obj.plt.lineSpec = p.Results.lineSpec;
            obj.plt.lineStyle = p.Results.lineStyle;
            obj.plt.color = p.Results.color;
            obj.plt.lineWidth = p.Results.lineWidth;
            obj.plt.marker = p.Results.marker;
            obj.plt.markerColor = p.Results.markerColor;
            obj.plt.name = p.Results.name;
            
            % Saving plot structure to the tile structure
            noOfPlots = length(obj.tiles{tileNo}.plots);
            obj.tiles{tileNo}.plots{noOfPlots+1} = obj.plt;
        end
        
        function obj = addSurf(obj, tileNo, srfArray, varargin)
            % Define default values
            defaultEdgeColor = 'auto';
            defaultLineStyle = '-';
            defaultLineWidth = 0.5;
            defaultFaceColor = 'auto';
            defaultFaceAlpha = 1;
            defaultName = 'none';
            
            % Define memberships sets
            lineStyleSet = {'-' '--' ':' '-.'};

            % Validation functions
            validTileNo = @(x) ge(x, 1) && isreal(x) && isnumeric(x) && eq(floor(x), x) && ge(obj.noOfTiles, x);
            validSrfArray = @(x) iscell(x) && ge(length(x), 3);
            validEdgeColor = @(x) validColorFunc(obj, x) || ismember(x, {'flat', 'interp'});
            validLineStyle = @(x) ismember(x, lineStyleSet);
            validLineWidth = @(x) ge(x, 0) && isreal(x) && isnumeric(x) && isscalar(x) && isfinite(x);
            validFaceColor = @(x) validColorFunc(obj, x) || ismember(x, {'flat', 'interp'});
            validFaceAlpha = @(x) ge(x, 0) && ge(0, x) && all(isreal(x)) && all(isscalar(x));
            validName = @(x) isstring(x) || ischar(x);
            
            % Validation function for each inner cell of srfArray -
            % assuming that the number length(srfArray) is greater than
            % three (potential gif frames).
            validGifFrame = @(x1, x2) iscell(x1) && iscell(x2) && eq(length(x1), length(x2)) && eq(length(x1), 3);
            
            % Input parser
            p = inputParser;
            p.CaseSensitive = true;

            % Adding arguments
            addRequired(p, 'tileNo', validTileNo);
            addRequired(p, 'srfArray', validSrfArray);
            addParameter(p, 'edgeColor', defaultEdgeColor, validEdgeColor);
            addParameter(p, 'lineStyle', defaultLineStyle, validLineStyle);
            addParameter(p, 'lineWidth', defaultLineWidth, validLineWidth);
            addParameter(p, 'faceColor', defaultFaceColor, validFaceColor);
            addParameter(p, 'faceAlpha', defaultFaceAlpha, validFaceAlpha);
            addParameter(p, 'name', defaultName, validName);

            parse(p, tileNo, srfArray, varargin{:});
            
            % Additional validation for tileNo
            if ~obj.definedTiles{p.Results.tileNo}
                error("This tile is not defined!");
            elseif obj.tiles{p.Results.tileNo}.dimensions == 2
                error("Tile must be three-dimensional!");
            end
            
            % Additional validation for the srfArray - check if the input
            % was specified for the gif animation
            if gt(length(p.Results.srfArray), 3)
                % It is potential gif input
                for i = 2:length(p.Results.srfArray)
                    % Check for validity of each cell according to
                    % validGifFrame function
                    if ~validGifFrame(p.Results.srfArray{i-1}, p.Results.srfArray{i})
                        error("Input srfArray is neither proper surf nor gif input.");
                    end
                end
                % Check if the number of frames is consistent with
                % other plots
                if gt(length(obj.tiles{tileNo}), 0) && ~eq(length(p.Results.srfArray), obj.tiles{tileNo}.plots{1}.noOfFrames)
                    error("Number of frames is too small/large for this tile.");
                else
                    % Set number of frames
                    obj.srf.noOfFrames = length(p.Results.srfArray);
                end
            else
                % Input is for a plot
                obj.srf.noOfFrames = 1;
            end

            % Extract variables from the parser
            tileNo = p.Results.tileNo;
            obj.srf.surfaceArray = p.Results.srfArray;
            obj.srf.edgeColor = p.Results.edgeColor;
            obj.srf.lineStyle = p.Results.lineStyle;
            obj.srf.lineWidth = p.Results.lineWidth;
            obj.srf.faceColor = p.Results.faceColor;
            obj.srf.faceAlpha = p.Results.faceAlpha;
            obj.srf.name = p.Results.name;
            
            % Saving surf structure to the tile structure
            noOfPlots = length(obj.tiles{tileNo}.plots);
            obj.tiles{tileNo}.plots{noOfPlots+1} = obj.srf;
        end
        
        function drawLayout(obj)
            % Looping through tiles
            for i = 1:length(obj.tiles)
                nexttile
                
                hold on
                obj.tiles{i}.handles = drawTile(obj, obj.tiles{i});
                hold off
                
                % Adding properties to each tile
                view(obj.tiles{i}.dimensions);
                
                if ~ismember(obj.tiles{i}.legend, {'none'})
                    % Exclude plot/surf 'none' value names from the legend 
                    names = [];
                    j = 1;
                    while j <= length(obj.tiles{i}.handles)
                        if isDefault(obj, obj.tiles{i}.plots{j}.name, 'none')
                            obj.tiles{i}.handles(j) = [];
                            j = j - 1;
                        else
                            names = [names, obj.tiles{i}.plots{j}.name];
                        end
                        j = j + 1;
                    end
                    
                    lgd = legend(obj.tiles{i}.handles);
                    lgd.String = names;
                    lgd.Location = obj.tiles{i}.legend;
                end
                
                if ~eq(obj.tiles{i}.title, 'none')
                    title(obj.tiles{i}.title)
                end
                
                xlabel(obj.tiles{i}.axesNames{1})
                ylabel(obj.tiles{i}.axesNames{2})
                if length(obj.tiles{i}.axesNames) == 3
                    zlabel(obj.tiles{i}.axesNames{3})
                end
                
                axis(obj.tiles{i}.size)
                
                grid on
            end
        end
    end
end