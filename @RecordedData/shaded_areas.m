function shaded_areas(x, y, error, varargin)
    % SHADED_AREAS Plots a line with shaded error areas.
    %   shaded_areas(x, y, error) plots the data in vectors x and y with
    %   shaded areas representing the error. The error can be a constant,
    %   a vector of the same length as y, or a matrix with two columns
    %   representing asymmetric errors.

    % Inputs:
    %   x - Vector of x data points
    %   y - Vector of y data points
    %   error - Error values (constant, column vector, or matrix with two columns for upper and lower bounds)
    %
    % Optional Name-Value Pair Inputs:
    %       'FaceColor' - Color of the shaded area (default: 'g')
    %       'FaceAlpha' - Transparency of the shaded area (default: 0.3)
    %       'LineWidth' - Width of the mean line (default: 3)
    %

    % Create an input parser
    p = inputParser;
    
    % Required parameters
    addRequired(p, 'x');
    addRequired(p, 'y');
    addRequired(p, 'error');
    
    % Optional parameter with default value
    addParameter(p, 'FaceColor', 'g');
    addParameter(p, 'FaceAlpha', 0.3);
    addParameter(p, 'LineWidth', 3);
    addParameter(p, 'DisplayName', []);
    
    % Parse inputs
    parse(p, x, y, error, varargin{:});
    
    % Extract values
    col = p.Results.FaceColor;
    alpha = p.Results.FaceAlpha;
    lw = p.Results.LineWidth;
    str2disp = p.Results.DisplayName;

    % Handle error
    if size(error,1) == length(y) & size(error,2) == 1 || numel(error) == 1 %error is a column vector the same length as y OR it is a constant value
        CIsup = y + error;
        CIinf = y - error;
    else size(error,1) == length(y) & size(error,2) == 2 %error is a matrix the same length as y with 2 columns
        CIsup = y + error(:,1);
        CIinf = y - error(:,2);
    end
    
    % Plot shaded area and line
    fill([x; flipud(x)], [CIsup; flipud(CIinf)], col, 'FaceAlpha', alpha, 'EdgeColor', col,'LineStyle', ':');
    hold on;
    if isempty(str2disp)
        plot(x, y, 'Color', col, 'LineWidth', lw);
    else
        plot(x, y, 'Color', col, 'LineWidth', lw, 'DisplayName', str2disp);
    end

end