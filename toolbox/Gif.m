% GIF - 
%
%   obj = Gif, constructor.
%           Output:
%       'obj': object of the class.
%
%   Limitations:
%       ...
%
%   Examples:
%       ...
%
%   Use:
%       ...
%
%   See also:
%       PARABOLOID, SPIRAL, CIRCLE, WAVEFUNCTION, QUANTUMN,
%       ENERGYAPPROXIMATION, WAVE, PLOT, DAVIDOVICRODS, CHANGENOTATIONTYPE
%
%   Patryk Jesionka, 2020
%
% Note: haven't tested for surfaces plotting

classdef Gif
    properties
        testLayout
        
        title
        axesNames
        size
        
        % Structure containing its pieces of information and plots
        % tile = struct('title', [], 'axesNames', [], 'size', [], 'dimensions', [], 'plots', {})
        % Structure containing plots data
        plt = struct('coordinatesArrays', [], 'lineStyle', [], 'lineColor', [], 'lineWidth', [], 'name', [], 'type', 'plot')
        srf = struct('surfaceArrays', [], 'edgeColor', [], 'lineStyle', [], 'faceColor', [], 'faceAlpha', [], 'name', [], 'type', 'surface')
        
        tile = []
        noOfFrames = 0
        noOfPlots = 0
    end
    methods
        function obj = defineTileGif(obj, title, axesNames, size)
            % Input validation using temporary Plot class object
            obj.testLayout = Plot;
            obj.testLayout = createLayout(obj.testLayout, 1, 1);
            obj.testLayout = defineTile(obj.testLayout, title, axesNames, size);
            
            % Assigning data to the class properties
            obj.title = title;
            obj.axesNames = axesNames;
            obj.size = size;
        end
        
        function obj = addPlotGif(obj, pltStruct, lStyle, lColor, lWidth, name)
            % Input validation
            obj.testLayout = addPlot(obj.testLayout, 1, pltStruct(1).coordinates, lStyle, lColor, lWidth, name);
            
            % Assigning data to the plot2d structure
            obj.plt.coordinatesArrays = pltStruct(1).coordinates;
            for i = 2:length(pltStruct)
                obj.plt.coordinatesArrays(i,:) = pltStruct(i).coordinates;
            end
            
            obj.plt.lineStyle = lStyle;
            obj.plt.lineColor = lColor;
            obj.plt.lineWidth = lWidth;
            obj.plt.name = name;
            
            if obj.noOfFrames == 0
                obj.noOfFrames = length(obj.plt.coordinatesArrays);
            else
                if ~(obj.noOfFrames == length(obj.plt.coordinatesArrays))
                    error("Invalid number of frames (number of coordinate arrays)");
                end
            end
            
            % Add plot to the tile
            obj.tile = [obj.tile obj.plt];
        end
        
        function obj = addSurfGif(obj, srfStruct, eColor, lStyle, fColor, fAlpha, name)
            % Input validation
            obj.testLayout = addSurf(obj, 1, srfStruct(1).coordinates, eColor, lStyle, fColor, fAlpha, name);
            
            % Assigning data to the plot3d structure
            for i = 1:length(srfStruct)
                obj.srf.surfaceArrays = [obj.srf.surfaceArrays srfStruct(i).coordinates];
            end
            obj.srf.edgeColor = eColor;
            obj.srf.lineStyle = lStyle;
            obj.srf.faceColor = fColor;
            obj.srf.faceAlpha = fAlpha;
            obj.srf.name = name;
            
            if obj.noOfFrames == 0
                obj.noOfFrames = length(obj.srf.surfaceArrays);
            else
                if ~(obj.noOfFrames == length(obj.srf.surfaceArrays))
                    error("Invalid number of frames (number of surface arrays)");
                end
            end
            
            % Ass surface to the tile
            obj.tile = [obj.tile obj.srf];
            
        end
        
        function drawGif(obj, filename)
            fig = figure;
            axis tight manual % this ensures that getframe() returns a consistent size
            if ~isstring(filename)
                error("filename must be a string!")
            end
            
            for i = 1:obj.noOfFrames

                layout = Plot;
                layout = createLayout(layout, 1, 1);
                layout = defineTile(layout, obj.title, obj.axesNames, obj.size);
                for j = 1:length(obj.tile)
                    if obj.tile(j).type == 'plot'
                        layout = addPlot(layout, 1, obj.tile(j).coordinatesArrays(i,:), obj.tile(j).lineStyle, obj.tile(j).lineColor, obj.tile(j).lineWidth, obj.tile(j).name);
                    else
                        layout = addSurf(layout, 1, obj.tile(j).surfacesArrays(i,:), obj.tile(j).endgeColor, obj.tile(j).lineStyle, obj.tile(j).faceColor, obj.tile(j).faceAlpha, obj.tile(j).name);
                    end
                end
                drawLayout(layout);

                % Capture the plot as an image
                frame = getframe(fig); 
                im = frame2im(frame); 
                [imind, cm] = rgb2ind(im, 256);

                % Write to the GIF File 
                if i == 1 
                    imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
                else 
                    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1); 
                end 
            end
        end
    end
end

% function gifEx
% %     % Create a circle in 2D (x-y plane) with radius 1 - 10
% %     circle1 = [];
% %     radius = 1:0.1:10;
% %     
% %     for i = radius
% %         circleHandle1 = Circle(radius(i));
% %         circle1 = [circle1 getCircle(circleHandle1)];
% %     end
% % 
% %     % % Create a circle in 2D (x-y plane) with radius 10 - 1
% %     circle2 = [];
% %     radius = 10:-0.1:1;
% %     
% %     for i = radius
% %         circleHandle2 = Circle(radius(i));
% %         circle2 = [circle2 getCircle(circleHandle2)];
% %     end
% % 
% %     title = "Resizing circles";
% %     axesNames = ["x" "y"];
% %     size = circle1(length(circle1)).size;
% % 
% %     pltStruct1 = circle1;
% %     lStyle1 = '-';
% %     lColor1 = 'k';
% %     lWidth1 = 0.5;
% %     name1 = "Circle #1";
% % 
% %     pltStruct2 = circle2;
% %     lStyle2 = '-';
% %     lColor2 = 'r';
% %     lWidth2 = 0.5;
% %     name2 = "Circle #2";
% % 
% %     fileName = "gifs/test.gif";
% % 
% %     gif = Gif;
% %     gif = defineTileGif(gif, title, axesNames, size);
% %     gif = addPlotGif(gif, pltStruct1, lStyle1, lColor1, lWidth1, name1);
% %     gif = addPlotGif(gif, pltStruct2, lStyle2, lColor2, lWidth2, name2);
% %     drawGif(gif, fileName);
%     
% %     %
% %     % Create a circle in 2D (x-y plane) with radius 0.1 - 10, centred at [0, 0]
% %     circle1 = [];
% %     changeVar1 = 1:0.1:10;
% %     for radius = changeVar1
% %         circleHandle1 = Circle(radius, [0 0], 'z', 360);
% %         circle1 = [circle1 getCircle(circleHandle1)];
% %     end
% % 
% %     % % Create a circle in 2D (x-y plane) with radius 10 - 0.1, centred at [0, 0]
% %     circle2 = [];
% %     changeVar2 = 10:-0.1:1;
% %     for radius = changeVar2
% %         circleHandle2 = Circle(radius, [0 0], 'z', 360);
% %         circle2 = [circle2 getCircle(circleHandle2)];
% %     end
% % 
% %     % input:
% %     title = "xyx";
% %     axesNames = ["x" "y"];
% %     size = graphs(length(graphs)).size;
% % 
% %     pltStruct1 = circle1;
% %     lStyle1 = '-';
% %     lColor1 = 'k';
% %     lWidth1 = 0.5;
% %     name1 = "Circle #1";
% % 
% %     pltStruct2 = circle2;
% %     lStyle2 = '-';
% %     lColor2 = 'r';
% %     lWidth2 = 0.5;
% %     name2 = "Circle #2";
% % 
% %     fileName = "gif/test.gif";
% % 
% %     % code:
% %     gif = Gif;
% %     gif = defineTileGif(gif, title, axesNames, size);
% %     gif = addPlotGif(gif, pltStruct1, lStyle1, lColor1, lWidth1, name1);
% %     gif = addPlotGif(gif, pltStruct2, lStyle2, lColor2, lWidth2, name2);
% %     drawGif(gif, fileName);
% end