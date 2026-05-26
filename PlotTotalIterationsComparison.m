%% Plot Total Iterations Comparison: Euler vs Leap-Frog
% This script reads total iterations data from Excel and creates bar charts
% comparing total iterations between Euler and Leap-Frog methods

function PlotTotalIterationsComparison()
    close all;
    
    % Specify the Excel file
    excelFile = "C:\Study\Y4S3\FYP\Total Iterations Extracted.xlsx";
    
    % Read data with original column names preserved
    data = readtable(excelFile, 'VariableNamingRule', 'preserve');
    
    % Extract CFL values from first column
    cfl_values = data{:, 1};
    
    % Define Reynolds numbers and data groups
    re_numbers = [100, 400, 1000];
    col_groups = {[1, 2, 3], [5, 6, 7], [9, 10, 11]};
    
    % Create figures for each Reynolds number
    for idx = 1:length(re_numbers)
        re = re_numbers(idx);
        cols = col_groups{idx};
        
        % Extract Euler and Leap-Frog data
        euler_iterations = data{:, cols(2)};
        lf_iterations = data{:, cols(3)};
        
        % Create figure
        fig = figure('Name', sprintf('Total Iterations Comparison - Re=%d', re), 'NumberTitle', 'off', "Theme", "light");
        fig.Position = [100 100 900 600];
        
        % Create bar chart
        x = categorical(string(cfl_values));
        y_data = [euler_iterations, lf_iterations];
        
        b = bar(x, y_data, 1.0);
        b(1).FaceColor = [0.2 0.4 0.8];      % Darker blue for Euler
        b(2).FaceColor = [0.8 0.2 0.2];      % Darker red for Leap-Frog
        
        % Set y-axis limits with padding
        max_val = max(y_data, [], 'all', 'omitnan');
        min_val = min(y_data, [], 'all', 'omitnan');
        margin = (max_val - min_val) * 0.2;
        bl = b.BaseLine;
        bl.BaseValue = min_val - margin;
        ylim([min_val - margin, max_val + margin]);
        
        xlabel('$\mathit{C}$ (Courant Number)', 'Interpreter','latex', 'FontSize',24);
        ylabel('Total Iterations', 'Interpreter', 'latex')
        legend('Euler', 'Leapfrog', 'Interpreter','latex', 'FontSize', 22, 'Location', 'northeast');
        grid on;
        set(gca, 'FontName', 'Cambria Math', 'FontSize', 22, 'XColor', 'k', 'YColor', 'k');
        
    end
    
end
