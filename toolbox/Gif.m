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
%       ENERGYAPPROXIMATION, WAVE, PLOT
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