function customizeSelectedArrowheads(h, G, edgesToModify, style)
    % Get the coordinates and sizes of nodes
    X = h.XData;
    Y = h.YData;
    nodeSizes = get(h, 'MarkerSize');  % Get the marker size (diameter in points)

    % Get axis limits for calculating scaling factors
    ax = ancestor(h, 'axes');
    xLimits = xlim(ax);
    yLimits = ylim(ax);
%     axPosition = get(ax, 'Position');  % Get axis position in figure
%     figPosition = get(get(ax,'Parent'), 'Position');  % Get figure position

%     % Convert marker size from points to data units
%     % Based on axis and figure size (accounting for aspect ratio)
%     scaleFactorX = diff(xLimits) / (axPosition(3) * figPosition(3) ); % 72 points per inch
%     scaleFactorY = diff(yLimits) / (axPosition(4) * figPosition(4) );

    currentunits = get(ax,'Units');
    set(ax, 'Units', 'Points');
    axpos = get(ax,'Position');
    set(ax, 'Units', currentunits);
    scaleFactorX = diff(xLimits)/axpos(3); % Calculate Marker width in points
    scaleFactorY = diff(yLimits)/axpos(4); % Calculate Marker width in points

    % Get the edge colors and line widths
    edgeColors = h.EdgeColor;
    lineWidths = h.LineWidth;
    
    hold(ax,'on'); % Ensure all customizations are kept in the plot

    % Loop through each edge to modify
    for i = 1:length(edgesToModify)
        edgeIndex = edgesToModify(i);
        [srcNode,tgtNode] = findedge(G,edgeIndex);
        
        % Get the positions of the source and target nodes
        x1 = X(srcNode);
        y1 = Y(srcNode);
        x2 = X(tgtNode);
        y2 = Y(tgtNode);
        
        % Get the radii of the source and target nodes in data units
        radiusSrc = (nodeSizes(srcNode) / 2) * sqrt(scaleFactorX^2 + scaleFactorY^2);
        radiusTgt = (nodeSizes(tgtNode) / 2) * sqrt(scaleFactorX^2 + scaleFactorY^2);

%         radiusSrc = sqrt(nodeSizes(srcNode)) * sqrt(scaleFactorX^2 + scaleFactorY^2) * 2;
%         radiusTgt = sqrt(nodeSizes(tgtNode)) * sqrt(scaleFactorX^2 + scaleFactorY^2) * 2;
        
        % Calculate the direction vector from the source to the target
        dx = x2 - x1;
        dy = y2 - y1;
        edgeLength = sqrt(dx^2 + dy^2);
        
        % Normalize the direction vector
        ux = dx / edgeLength;
        uy = dy / edgeLength;
        
        % Adjust the start and end points based on node radii
        x1Adjusted = x1 + ux * radiusSrc;
        y1Adjusted = y1 + uy * radiusSrc;
        x2Adjusted = x2 - ux * radiusTgt;
        y2Adjusted = y2 - uy * radiusTgt;

        % Determine the color for the current edge
        if size(edgeColors, 1) == 1
            % Single color for all edges
            edgeColor = edgeColors;
        else
            % Different colors for each edge
            edgeColor = edgeColors(edgeIndex, :);
        end

        if length(lineWidths) == 1
            % Single line width for all edges
            lineWidth = lineWidths;
        else
            % Different line widths for each edge
            lineWidth = lineWidths(edgeIndex);
        end

        % Set the arrowhead or perpendicular line based on style
        switch style
            case 'arrow'
                % Standard arrowhead
                quiver(ax,x1Adjusted, y1Adjusted, x2Adjusted - x1Adjusted, y2Adjusted - y1Adjusted, 0, ...
                    'MaxHeadSize', 1, 'Color', edgeColor, 'LineWidth', lineWidth, 'ShowArrowHead', true, 'Autoscale', 'off');
                
            case 'perpendicular'
                % Draw a perpendicular line at the target node
                % Draw the main edge line
                line(ax,[x1Adjusted, x2Adjusted], [y1Adjusted, y2Adjusted], 'Color', edgeColor, 'LineWidth', lineWidth);
                
                % Calculate the length of the perpendicular line
                len = 3*lineWidth* sqrt(scaleFactorX^2 + scaleFactorY^2);%0.1 * edgeLength;%3*lineWidth* sqrt(scaleFactorX^2 + scaleFactorY^2);%0.05 * edgeLength;
                
                % Perpendicular line direction
                perpAngle = atan2(dy, dx) + pi/2;
                px = len * cos(perpAngle);
                py = len * sin(perpAngle);
                
                % Draw the perpendicular line
                line(ax,[x2Adjusted - px, x2Adjusted + px], [y2Adjusted - py, y2Adjusted + py], 'Color', edgeColor, 'LineWidth', lineWidth);
                
            otherwise % no head
                % Draw a simple line (no arrowhead)
                line(ax,[x1Adjusted, x2Adjusted], [y1Adjusted, y2Adjusted], 'Color', edgeColor, 'LineWidth', lineWidth);
        end
    end

    hold(ax, 'off'); % Release the plot hold
end
