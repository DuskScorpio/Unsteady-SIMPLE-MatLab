%% Plot Euler Simulation Results vs Ghia Reference Data
% This script reads CSV files from a simulation results directory and compares
% them with Ghia et al. (1982) benchmark values

function PlotEulerResults()
    close all;
    
    % Get script directory for Ghia data
    scriptDir = fileparts(mfilename('fullpath'));
    ghiaFolder = fullfile(scriptDir, 'Ghia Centerline');
    
    % Read Ghia et al. (1982) benchmark values
    u_table = readtable(fullfile(ghiaFolder, 'Ghia-u-centerline.csv'), ReadVariableNames=true, VariableNamingRule="preserve");
    v_table = readtable(fullfile(ghiaFolder, 'Ghia-v-centerline.csv'), ReadVariableNames=true, VariableNamingRule="preserve");
    
    % Define Reynolds numbers and corresponding directories
    re_values = [100, 400, 1000];
    euler_dirs = {
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Euler\Results rev1\Re=100 CFL=1e-01 15-Apr_21-12-53";
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Euler\Results rev1\Re=400 CFL=1e-01 15-Apr_21-20-32";
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Euler\Results rev1\Re=1000 CFL=1e-01 15-Apr_21-49-10"
    };
    leapfrog_dirs = {
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Leap-Frog\Results rev1\Re=100 CFL=1e-01 15-Apr_22-28-17";
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Leap-Frog\Results rev1\Re=400 CFL=1e-01 15-Apr_22-43-18";
        "C:\Study\Y4S3\FYP\MatLab\Results_LidDrivenCavity_Leap-Frog\Results rev1\Re=1000 CFL=1e-01 15-Apr_23-10-37"
    };
    
    % Plot for each Reynolds number
    for idx = 1:length(re_values)
        re = re_values(idx);
        eulerDir = euler_dirs{idx};
        leapfrogDir = leapfrog_dirs{idx};
        
        % Build file names
        u_file_euler = sprintf('Re=%d CFL=1e-01 u_Center_Line.csv', re);
        v_file_euler = sprintf('Re=%d CFL=1e-01 v_Center_Line.csv', re);
        u_file_lf = sprintf('Re=%d CFL=1e-01 u_Center_Line.csv', re);
        v_file_lf = sprintf('Re=%d CFL=1e-01 v_Center_Line.csv', re);
        
        % Read Euler simulation results from CSV files
        u_data_euler = readmatrix(fullfile(eulerDir, u_file_euler));
        v_data_euler = readmatrix(fullfile(eulerDir, v_file_euler));
        
        % Read Leap-Frog simulation results from CSV files
        u_data_lf = readmatrix(fullfile(leapfrogDir, u_file_lf));
        v_data_lf = readmatrix(fullfile(leapfrogDir, v_file_lf));
        
        % Reconstruct grid coordinates based on actual data size
        num_u = length(u_data_euler);
        num_v = length(v_data_euler);
        
        % Create coordinate arrays matching data size
        y_center = linspace(1, 0, num_u)';
        u_centerline = u_data_euler;
        u_centerline_lf = u_data_lf;
        
        x_center = linspace(0, 1, num_v)';
        v_centerline = v_data_euler;
        v_centerline_lf = v_data_lf;
        
        % Extract Ghia data for current Re
        ghia_y = u_table.y;
        ghia_u = u_table.(sprintf("%d", re));
        
        ghia_x = v_table.x;
        ghia_v = v_table.(sprintf("%d", re));
        
        % Create figure with dual y-axes
        figure("Theme", "light");
        
        % Left y-axis for u velocity
        yyaxis left
        h1 = plot(y_center, u_centerline, 'k', 'LineWidth', 1.5); 
        hold on;
        h2 = plot(y_center, u_centerline_lf, 'g--', 'LineWidth', 1.5); 
        h3 = plot(ghia_y, ghia_u, 'rs', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); 
        ylabel('$u$', 'Interpreter', 'latex');
        ax1 = gca;
        ax1.FontName = 'Cambria Math';
        ax1.YAxis(1).Color = 'k';
        ax1.XAxis.Color = 'k';
        
        % Right y-axis for v velocity
        yyaxis right
        plot(x_center, v_centerline, 'k', 'LineWidth', 1.5); 
        hold on;
        plot(x_center, v_centerline_lf, 'g--', 'LineWidth', 1.5); 
        h4 = plot(ghia_x, ghia_v, 'b^', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        ylabel('$v$', 'Interpreter', 'latex');
        ax2 = gca;
        ax2.FontName = 'Cambria Math';
        ax2.YAxis(2).Color = 'k';
        ax2.XAxis.Color = 'k';
        
        xlabel('$X(Y)$', 'Interpreter', 'latex');
        set(gca, 'FontName', 'Cambria Math', 'XColor', 'k', 'YColor', 'k');
        
        % Only show legend for Re=1000
        if idx == 3
            h_legend = legend([h1, h3, h2, h4], 'Euler', 'Ghia $u$', 'Leapfrog','Ghia $v$', 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2, 'edgecolor', 'none', 'Interpreter', 'latex');
            h_legend.FontName = 'Cambria Math';
            h_legend.FontSize = 11;
        end
        axis square;
        grid off;
        
        % Calculate and display RMSE for reference
        % Interpolate simulation data to Ghia grid points
        u_centerline_interp = interp1(y_center, u_centerline, ghia_y, 'linear');
        v_centerline_interp = interp1(x_center, v_centerline, ghia_x, 'linear');
        u_centerline_interp_lf = interp1(y_center, u_centerline_lf, ghia_y, 'linear');
        v_centerline_interp_lf = interp1(x_center, v_centerline_lf, ghia_x, 'linear');
        
        u_rmse = sqrt(mean((u_centerline_interp - ghia_u).^2, "omitnan"));
        v_rmse = sqrt(mean((v_centerline_interp - ghia_v).^2, "omitnan"));
        u_rmse_lf = sqrt(mean((u_centerline_interp_lf - ghia_u).^2, "omitnan"));
        v_rmse_lf = sqrt(mean((v_centerline_interp_lf - ghia_v).^2, "omitnan"));
        
        fprintf('Re=%d - RMSE for u-velocity (Euler): %.6f\n', re, u_rmse);
        fprintf('Re=%d - RMSE for v-velocity (Euler): %.6f\n', re, v_rmse);
        fprintf('Re=%d - RMSE for u-velocity (Leap-Frog): %.6f\n', re, u_rmse_lf);
        fprintf('Re=%d - RMSE for v-velocity (Leap-Frog): %.6f\n', re, v_rmse_lf);
    end
end
