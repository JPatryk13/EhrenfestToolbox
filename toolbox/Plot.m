%
%
%      INSTRUCTION
% Call in the order:
% 1. createLayout
% 2. defineTile (as many times as you need)
% 3. addPlot, addSurf (as many times as you need)
% 4. drawLayout 
% 
%

classdef Plot
    properties
        % Tiled layout handle 
        layout
        % Number of tiles available - used to validate if user specified
        % existing tile_no
        noOfTiles {mustBePositive} = 1
        % Contains information whether the certain tile was defined or not
        definedTiles = {}
        
        % Array containing structures for every tile
        tiles = {}
        % Structure containing its pieces of information and plots
        % tile = struct('title', [], 'axesNames', [], 'size', [], 'dimensions', [], 'plots', {})
        % Structure containing plots data
        plt = struct('coordinatesArray', [], 'lineStyle', [], 'lineColor', [], 'lineWidth', [], 'name', [], 'type', 'plot')
        srf = struct('syrfaceArray', [], 'edgeColor', [], 'lineStyle', [], 'faceColor', [], 'faceAlpha', [], 'name', [], 'type', 'surface')
        
        % Plot properties - used to validate user input and set default
        % values.
        lineStyle {mustBeMember(lineStyle, {'-', '--', ':', '-.'})} = '-'
        lineColor {mustBeMember(lineColor, {'r', 'g', 'b', 'w', 'c', 'm', 'y', 'k'})} = 'k'
        lineWidth {mustBePositive} = 0.5
        
        % Surface properties - used to validate user input
        edgeColor {mustBeMember(edgeColor, {'none', 'flat'})} = 'none'
        faceColor = 'interp'
        faceAlpha {mustBeGreaterThanOrEqual(faceAlpha, 0), mustBeLessThanOrEqual(faceAlpha, 1)} = 1
    end
    methods (Access = private)
        function drawPlot(~, object, dimensions)
            % Creating plot object (first assigning coordinates for the
            % sake of explicity)
            x = object.coordinatesArray{1};
            y = object.coordinatesArray{2};
            if dimensions == 3
                z = object.coordinatesArray{3};
                p = plot3(x, y, z);
            else
                p = plot(x, y);
            end
            
            % Setting properties
            p.LineStyle = object.lineStyle;
            p.Color = object.lineColor;
            p.LineWidth = object.lineWidth;
            p.DisplayName = object.name;
        end
        
        function drawSurf(~, object)
            % Creating surface object (first assigning coordinates for the
            % sake of explicity)
            X = object.surfaceArray{1};
            Y = object.surfaceArray{2};
            Z = object.surfaceArray{3};
            s = surf(X, Y, Z);
            % Setting properties
            s.EdgeColor = object.edgeColor;
            s.LineStyle = object.lineStyle;
            s.FaceColor = object.faceColor;
            s.FaceAlpha = object.faceAlpha;
            s.DisplayName = object.name;
        end
        
        function drawTile(obj, tile)
            % Looping through plots
            for j = 1:length(tile.plots)
                % If it is a 2d plot
                if tile.dimensions == 2
                    drawPlot(obj, tile.plots{j}, 2);
                elseif tile.dimensions == 3 % if it is three-dimensional
                    if tile.plots{j}.type == "plot"
                        drawPlot(obj, tile.plots{j}, 3);
                    elseif tile.plots{j}.type == "surface"
                        drawSurf(obj, tile.plots{j});
                    end
                end
            end
        end
    end
    methods (Access = public)
        function obj = createLayout(obj, m, n)
            % Creates m by n layout; then, specifies number of tiles
            obj.layout = tiledlayout(m, n);
            obj.noOfTiles = m*n;
            
            % Iterates through the number of tiles and adds false values 
            % since tiles are not defined yet
            for i = 1:obj.noOfTiles
                obj.definedTiles{i} = false;
            end
        end
        
        function obj = defineTile(obj, title, axesNames, size)
            % title = "title"
            % axesNames = {'x' 'y' 'z'}
            % sizes = [-x x -y y -z z]
            
            % Input validation
            if ~isstring(title)
                error("Title must be string type!");
            elseif ~(length(axesNames) == 0.5*length(size))
                error("Specified size and axesNames arrays does not cover the same dimensionality!");
            elseif ~(2 >= length(axesNames) <= 3)
                error("axesNames array must have either 2 or 3 dimensions!");
            elseif ~(length(size) == 4 || length(size) == 6)
                error("size array must have either 4 or 6 dimensions!");
            end
            % Check if array content has a valid type
            for i = 1:length(axesNames)
                if ~(isstring(axesNames{i}) || ischar(axesNames{i}))
                    error("All axes names must be string type!")
                end
            end
            for i = 1:length(size)
                if ~isnumeric(size(i))
                    error("Size must be an array with integers!")
                end
            end
            % Check if there is a space for another tile
            if ~(length(obj.tiles)+1 <= obj.noOfTiles)
                error("There is no space for another tile in the layout!");
            end
            
            % Assigning data to the tile structure
            tile.title = title;
            tile.axesNames = axesNames;
            tile.size = size;
            tile.dimensions = length(axesNames);
            tile.plots = {};
            
            % Adds structure to the tiles array
            obj.tiles{length(obj.tiles)+1} = tile;
            
            % Set the first false flag in definedTiles to true
            for i = 1:obj.noOfTiles
                if ~obj.definedTiles{i}
                    obj.definedTiles{i} = true;
                end
            end
        end
        
        function obj = addPlot(obj, tileNo, pltArray, lStyle, lColor, lWidth, name)
            % Input validation
            if tileNo <= 0 && tileNo >= obj.noOfTiles
                error("Invalid tile number specified!");
            elseif ~obj.definedTiles{tileNo}
                error("That tile is not defined!");
            elseif ~(length(pltArray) == 2 || length(pltArray) == 3)
                error("pltArray must be two or three-dimensional!");
            elseif ~(obj.tiles{tileNo}.dimensions == length(pltArray))
                error("Tile has different dimensions than the plot you want to add!");
            elseif ~isstring(name)
                error("Name of the plot must be a string type!");
            end
            obj.lineStyle = lStyle;
            obj.lineColor = lColor;
            obj.lineWidth = lWidth;
            
            % Assigning data to the plot2d structure
            obj.plt.coordinatesArray = pltArray;
            obj.plt.lineStyle = obj.lineStyle;
            obj.plt.lineColor = obj.lineColor;
            obj.plt.lineWidth = obj.lineWidth;
            obj.plt.name = name;
            
            % Saving plot structure to the tile structure
            noOfPlots = length(obj.tiles{tileNo}.plots);
            obj.tiles{tileNo}.plots{noOfPlots+1} = obj.plt;
        end
        
        function obj = addSurf(obj, tileNo, srfArray, eColor, lStyle, fColor, fAlpha, name)
            % Input Validation
            if tileNo <= 0 && tileNo >= obj.noOfTiles
                error("Invalid tile number specified!");
            elseif ~obj.definedTiles{tileNo}
                error("That tile is not defined!");
            elseif ~(length(srfArray) == 3)
                error("pltArray must be three-dimensional!");
            elseif obj.tiles{tileNo}.dimensions == 2
                error("Tile has two dimensions while the plot you want to add is three-dimensional!");
            elseif ~isstring(name)
                error("Name of the plot must be a string type!");
            end
            obj.edgeColor = eColor;
            obj.lineStyle = lStyle;
            obj.faceColor = fColor;
            obj.faceAlpha = fAlpha;
            
            % Assigning data to the plot3d structure
            obj.srf.surfaceArray = srfArray;
            obj.srf.edgeColor = obj.edgeColor;
            obj.srf.lineStyle = obj.lineStyle;
            obj.srf.faceColor = obj.faceColor;
            obj.srf.faceAlpha = obj.faceAlpha;
            obj.srf.name = name;
            
            % Saving surf structure to the tile structure
            noOfPlots = length(obj.tiles{tileNo}.plots);
            obj.tiles{tileNo}.plots{noOfPlots+1} = obj.srf;
        end
        
        function drawLayout(obj)
            % Looping through tiles
            for i = 1:length(obj.tiles)
                nexttile
                hold on
                drawTile(obj, obj.tiles{i});
                hold off
                
                % Adding properties to each tile
                legend
                view(obj.tiles{i}.dimensions)
                title(obj.tiles{i}.title)
                xlabel(obj.tiles{i}.axesNames{1})
                ylabel(obj.tiles{i}.axesNames{2})
                axis(obj.tiles{i}.size)
                if length([]) == 3
                    zlabel([obj.tiles{i}.axesNames{3}])
                end
                grid on
            end
        end
        
%         function gifLayout()
%             n_laps = 0:0.05:4.95;
%             nImages = length(n_laps);
% 
%             n_rmin = 4.95:-0.05:0;
% 
%             fig = figure;
%             for idx = 1:nImages
%                 % hold on
%                 % createParaboloid(4, 1)
%                 createCircle(0, 0, 0, 5, n_rmin(idx), n_laps(idx), [1 0], 5, true, true, true)
%                 % hold off
%                 drawnow
%                 frame = getframe(fig);
%                 im{idx} = frame2im(frame);
%             end
% 
%             filename = 'm3d_01.gif'; % Specify the output file name
%             for idx = 1:nImages
%                 [A,map] = rgb2ind(im{idx},256);
%                 if idx == 1
%                     imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
%                 else
%                     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%                 end
%             end
%             close;
%         end
    end
end