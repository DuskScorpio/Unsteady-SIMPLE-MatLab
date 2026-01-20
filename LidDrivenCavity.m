clearvars -except Re
close all
clc

%% Setup
% Primary Parameters
length = 1; % Length along the positive x-directon of the flow domain 
height = 1; % Length along the positive y-direction of the flow domain
numCellsY = 51; % Number of cells along the y-direction
CFL = 0.5; % Courant number
U_lid = 1; % Lid velocity = 1m/s
alpha_v = 1; % v-velocity relaxation factor
alpha_u = 1; % u-velocity relaxation factor
alpha_p = 0.3; % pressure relaxation factor

% Derived Parameters
dy = height / numCellsY; % Cell size along the y-direction
dx = dy; % Cell size along the x-direction
numCellsX = length / dx; % Number of cells along the x-direction
dt = CFL * (dx / U_lid); % Time step size

% Criteria
timeEnd = 100; % Termination time
maxSteps = round(timeEnd / dt); % Max number of time steps
maxIterations = 200; % Max number of SIMPLE iterations
maxResidual = 1e-8; % SIMPLE loop residual tolerance
steadyTolerance = 1e-6; % steady state tolerance

if ~exist('Re', 'var')
    Re = 100; % Set Re if script ran by itself
end

% Ensure uniform grid size
if mod(numCellsX,1) ~= 0
    error('numCellsX is not an integer. Adjust numCellsY or dimensions.');
end

% Fluid properties
MU = 1e-3; % Viscosity of the fluid
RHO = (MU * Re) / (U_lid * length); % Density of the fluid

%% Variables
u_old = zeros(numCellsY + 2, numCellsX + 3); % u velocity from previous time step
u_star = zeros(numCellsY + 2, numCellsX + 3); % u velocity prediction
u_new = zeros(numCellsY + 2, numCellsX + 3); % Corrected u velocity

v_old = zeros(numCellsY + 3, numCellsX + 2); % v velocity from previous time step
v_star = zeros(numCellsY + 3, numCellsX + 2); % v velocity prediction
v_new = zeros(numCellsY + 3, numCellsX + 2); % Corrected v velocity

p_prime = zeros(numCellsY + 2, numCellsX + 2); % Pressure correction
p_old = zeros(numCellsY + 2, numCellsX + 2); % Pressure from previous time step
p = zeros(numCellsY + 2, numCellsX + 2); % Pressure
p_new = zeros(numCellsY + 2, numCellsX + 2); % Corrected pressure

b = zeros(numCellsY + 2, numCellsX + 2); % Source term (for pressure correction)
b_new = zeros(numCellsY + 2, numCellsX + 2); % Source term (after correction)

%% Initial Conditions
u_old(1, 2:numCellsX + 2) = 2*U_lid; % u velocity at lid = 1m/s
u = u_old; % u velocity from previous iteration
v = v_old; % v velocity from previous iteration

steadyReached = false;
n = 0; % time step count
iterations = 0;
totalIterations = 0;
residual = 1;
maxDiff = 1;

%% Build the sparse matrix A for the pressure correction equation
Nx = numCellsX;
Ny = numCellsY;
N  = Nx * Ny;

a_EW = dt / (RHO * dx * dx);
a_NS = dt / (RHO * dy * dy);
a_P  = 2*a_EW + 2*a_NS;

A = sparse(N, N);

for i = 1:Ny
    for j = 1:Nx
        k = (i - 1) * Nx + j; % Row major index

        A(k,k) = a_P; % Center

        if j > 1
            A(k, k-1) = -a_EW; % West neighbour
        else
            A(k,k) = A(k,k) - a_EW; % West boundary
        end

        if j < Nx
            A(k, k+1) = -a_EW; % East neighbour
        else
            A(k,k) = A(k,k) - a_EW; % East boundary
        end

        if i < Ny
            A(k, k+Nx) = -a_NS; % South neighbour
        else
            A(k,k) = A(k,k) - a_NS; % South boundary
        end

        if i > 1
            A(k, k-Nx) = -a_NS; % North neighbour
        else
            A(k,k) = A(k,k) - a_NS; % North boundary
        end
    end
end

% Fix singularity
A(1,:) = 0;
A(1,1) = 1;

%% Residual Plot Setup
set(0,'DefaultFigureWindowStyle','docked')
figure;
hResidual = animatedline('Color', 'red', 'LineWidth', 1.5, 'DisplayName', 'Mass Residual (Inner Loop)');
hUDiff = animatedline('Color', 'green', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'u Diff (Outer Loop)');
hVDiff = animatedline('Color', 'blue', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'v Diff (Outer Loop)');
hPDiff = animatedline('Color', 'magenta', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'p Diff (Outer Loop)');
set(gca, 'YScale', 'log');
set(gca, 'YScale', 'log');
grid on;
axis square;
xlabel('Total SIMPLE Iterations');
ylabel('RMS Residual');
title('Solver Convergence');
legend('show', 'Location', 'southwest');

%% Start timer
tic

%% Time marching
while ~steadyReached && n < maxSteps

    % Print current time, iteration, and residual
    currentTime = n * dt; % if inside the outer time loop
    fprintf('Time = %.3f s, Iteration = %d, Residual = %.6e, Tolerance = %.6e\n', currentTime, iterations, residual, maxDiff);

    n = n + 1;
    iterations = 0;
    residual = 1;
    
    %% SIMPLE loop
    while residual > maxResidual && iterations < maxIterations
    
        iterations = iterations + 1;
    
        %% x-momentum equation - Interior
        for i = 2:numCellsY + 1
            for j = 3:numCellsX + 1
    
                u_P = u_old(i, j);
                u_E = u_old(i, j + 1);
                u_EE = u_old(i, j + 2);
                u_W = u_old(i, j - 1);
                u_WW = u_old(i, j - 2);
                u_N = u_old(i - 1, j);
                u_S = u_old(i + 1, j);
    
                u_e_interp = 0.5 * (u_P + u_E);
                u_w_interp = 0.5 * (u_P + u_W);
                u_n_interp = 0.5 * (u_P + u_N);
                u_s_interp = 0.5 * (u_P + u_S);
                v_n_interp = 0.5 * (v_old(i, j) + v_old(i, j - 1));
                v_s_interp = 0.5 * (v_old(i + 1, j) + v_old(i + 1, j - 1));
    
                p_w = p(i, j - 1);
                p_e = p(i, j);
    
                % Second-order upwind for u_e
                if u_e_interp < 0
                    u_e = 1.5*u_E - 0.5*u_EE;
                elseif u_e_interp > 0
                    u_e = 1.5*u_P - 0.5*u_W;
                else
                    u_e = 0;
                end
    
                % Second-order upwind for u_w
                if u_w_interp < 0
                    u_w = 1.5*u_P - 0.5*u_E;
                elseif u_w_interp > 0
                    u_w = 1.5*u_W - 0.5*u_WW;
                else
                    u_w = 0;
                end
    
                convectionTerm = (u_w_interp * u_w - u_e_interp * u_e) / dx + (v_s_interp * u_s_interp - v_n_interp * u_n_interp) / dy; % LUDS
                diffusionTerm = (MU / RHO) * ((u_E - 2*u_P + u_W) / (dx * dx) + (u_N - 2*u_P + u_S) / (dy * dy));
                pressureTerm = (p_w - p_e) / (RHO * dx);
    
                u_star(i, j) = 0.25*(u_e + u_w + u_n_interp + u_s_interp) + dt * (convectionTerm + diffusionTerm + pressureTerm);
                
            end
        end
    
        %% x-momentum equation - Boundary
        % Left wall (u = 0)
        u_star(2:numCellsY + 1, 2) = 0; 
    
        % Left ghost cell
        u_star(2:numCellsY + 1, 1) = - u_star(2:numCellsY + 1, 3);
    
        % Right wall (u = 0)
        u_star(2:numCellsY + 1, numCellsX + 2) = 0;
    
        % Right ghost cell
        u_star(2:numCellsY + 1, numCellsX + 3) = - u_star(2:numCellsY + 1, numCellsX + 1);
    
        % Top wall (u = 1)
        u_star(1, 2:numCellsX + 2) = 2*U_lid - u_star(2, 2:numCellsX + 2);
    
        % Bottom wall (u = 0)
        u_star(numCellsY + 2, 2:numCellsX + 2) = -(u_star(numCellsY + 1, 2:numCellsX + 2));
    
        %% y-momentum eq. - Interior
        for i = 3:numCellsY + 1
            for j = 2:numCellsX + 1
                         
                v_P = v_old(i, j);
                v_E = v_old(i, j + 1);
                v_W = v_old(i, j - 1);
                v_N = v_old(i - 1, j);
                v_NN = v_old(i - 2, j);
                v_S = v_old(i + 1, j);
                v_SS = v_old(i + 2, j);
    
                u_e_interp = 0.5 * (u_old(i, j + 1) + u_old(i - 1, j + 1));
                u_w_interp = 0.5 * (u_old(i, j) + u_old(i - 1, j));
                v_e_interp = 0.5 * (v_P + v_E);
                v_w_interp = 0.5 * (v_P + v_W);
                v_n_interp = 0.5 * (v_P + v_N);
                v_s_interp = 0.5 * (v_P + v_S);
    
                p_n = p(i - 1, j);
                p_s = p(i, j);
    
                % Second-order upwind for v_n
                if v_n_interp < 0
                    v_n = 1.5*v_N - 0.5*v_NN;
                elseif v_n_interp > 0
                    v_n = 1.5*v_P - 0.5*v_S;
                else
                    v_n = 0;
                end
    
                % Second-order upwind for v_s
                if v_s_interp < 0
                    v_s = 1.5*v_P - 0.5*v_N;
                elseif v_s_interp > 0
                    v_s = 1.5*v_S - 0.5*v_SS;
                else
                    v_s = 0;
                end
    
                convectionTerm = (u_w_interp * v_w_interp - u_e_interp * v_e_interp) / dx + (v_s_interp * v_s - v_n_interp * v_n) / dy; % LUDS
                diffusionTerm = (MU / RHO) * ((v_E - 2*v_P + v_W) / (dx * dx) + (v_N - 2*v_P + v_S) / (dy * dy));
                pressureTerm = (p_s - p_n) / (RHO * dy);
                v_star(i, j) = 0.25*(v_e_interp + v_w_interp + v_n + v_s) + dt * (convectionTerm + diffusionTerm + pressureTerm);
    
            end
        end
    
        %% v-momentum eq. - Boundary
        % Top wall (v = 0)
        v_star(2, 2:numCellsX + 1) = 0;
    
        % Top ghost cell
        v_star(1, 2:numCellsX + 1) = - v_star(3, 2:numCellsX + 1);
    
        % Bottom wall (v = 0)
        v_star(numCellsY + 2, 2:numCellsX + 1) = 0;
    
        % Bottom ghost cell
        v_star(numCellsY + 3, 2:numCellsX + 1) = - v_star(numCellsY + 1, 2:numCellsX + 1);
    
        % Left wall (v = 0)
        v_star(2:numCellsY + 2, 1) = - v_star(2:numCellsY + 2, 2);
    
        % Right wall (v = 0)
        v_star(2:numCellsY + 2, numCellsX + 2) = - v_star(2:numCellsY + 2, numCellsX + 1);

        %% Pressure correction - Interior
        for i = 2:numCellsY + 1
            for j = 2:numCellsX + 1
                b(i, j) = (u_star(i, j) - u_star(i, j + 1)) / dx + (v_star(i + 1, j) - v_star(i, j)) / dy;
            end
        end

        % Build RHS
        b_vec = zeros(N,1);

        for i = 1:Ny
            for j = 1:Nx
                k = (i - 1) * Nx + j;
                b_vec(k) = b(i + 1, j + 1);
            end
        end

        % Fix singularity
        b_vec(1) = 0;

        % Solve
        p_vec = A \ b_vec;

        % Map back to staggered grid p_prime(i,j)
        p_prime(:) = 0;   % reset
        for i = 1:Ny
            for j = 1:Nx
                k = (i - 1) * Nx + j;
                p_prime(i + 1, j + 1) = p_vec(k);
            end
        end
    
        %% Correct pressures
        for j = 2:numCellsX + 1
            for i = 2:numCellsY+1
    
                p_new(i, j) = p(i, j) + alpha_p * p_prime(i, j);
    
            end
        end
    
        %% Pressure field - Boundary
        % Left wall
        p_new(2:numCellsY + 1, 1) = p_new(2:numCellsY + 1, 2);
    
        % Right wall
        p_new(2:numCellsY + 1, numCellsX + 2) = p_new(2:numCellsY + 1, numCellsX + 1);
    
        % Top wall
        p_new(1, 2:numCellsX + 1) = p_new(2, 2:numCellsX + 1);
    
        % Bottom wall
        p_new(numCellsY + 2, 2:numCellsX + 1) = p_new(numCellsY + 1, 2:numCellsX + 1);
    
        %% Correct u velocities - Interior
        for i = 2:numCellsY + 1
            for j = 3:numCellsX + 1
           
            u_new(i, j) = u_star(i, j) + alpha_u * dt / (RHO * dx) * (p_prime(i, j - 1) - p_prime(i, j));
    
            end
        end
    
        %% Correct u velocities - Boundary
        % Left wall (u = 0)
        u_new(2:numCellsY + 1, 2) = 0; 
    
        % Left ghost cell
        u_new(2:numCellsY + 1, 1) = - u_new(2:numCellsY + 1, 3);
    
        % Right wall (u = 0)
        u_new(2:numCellsY + 1, numCellsX + 2) = 0;
    
        % Right ghost cell
        u_new(2:numCellsY + 1, numCellsX + 3) = - u_new(2:numCellsY + 1, numCellsX + 1);
    
        % Top wall (u = 1)
        u_new(1, 2:numCellsX + 2) = 2*U_lid - u_new(2, 2:numCellsX + 2);
    
        % Bottom wall (u = 0)
        u_new(numCellsY + 2, 2:numCellsX + 2) = -(u_new(numCellsY + 1, 2:numCellsX + 2));
    
        %% Correct v velocities - Interior
        for i = 3:numCellsY + 1
            for j = 2:numCellsX + 1
    
                v_new(i, j) = v_star(i, j) + alpha_v * dt / (RHO * dy) * (p_prime(i + 1, j) - p_prime(i, j));
    
            end
        end
    
        %% Correct v velocities - Boundary
        % Top wall (v = 0)
        v_new(2, 2:numCellsX + 1) = 0;
    
        % Top ghost cell
        v_new(1, 2:numCellsX + 1) = - v_new(3, 2:numCellsX + 1);
    
        % Bottom wall (v = 0)
        v_new(numCellsY + 2, 2:numCellsX + 1) = 0;
    
        % Bottom ghost cell
        v_new(numCellsY + 3, 2:numCellsX + 1) = - v_new(numCellsY + 1, 2:numCellsX + 1);
    
        % Left wall (v = 0)
        v_new(2:numCellsY + 2, 1) = - v_new(2:numCellsY + 2, 2);
    
        % Right wall (v = 0)
        v_new(2:numCellsY + 2, numCellsX + 2) = - v_new(2:numCellsY + 2, numCellsX + 1);
    
        %% Compute residual
        for i = 2:numCellsY + 1
            for j = 2:numCellsX + 1
    
                b_new(i, j) = dy * (u_new(i, j) - u_new(i, j + 1)) + dx * (v_new(i + 1, j) - v_new(i, j));
    
            end
        end
    
        residual = sqrt(mean(b_new.^2, 'all')); % compute residual as RMS of the source term

        % Update residual graph
        addpoints(hResidual, totalIterations + iterations, residual);
        
        % Update the plot every 5 iterations to keep the solver fast
        if mod(totalIterations + iterations, 100) == 0
            xlim([max(1, totalIterations + iterations - 500), totalIterations + iterations + 50]); % Show latest 500 iterations
            drawnow limitrate
        end
    
        %% End of iteration
        u = u_new;
        v = v_new;
        p = p_new;
    
    end

    %% End of time step
    % Verify steady state
    u_diff = sqrt(mean((u - u_old).^2, 'all'));
    v_diff = sqrt(mean((v - v_old).^2, 'all'));
    p_diff = sqrt(mean((p - p_old).^2, 'all'));
    maxDiff = max(u_diff, v_diff);

    addpoints(hUDiff, totalIterations, u_diff);
    addpoints(hVDiff, totalIterations, v_diff);
    addpoints(hPDiff, totalIterations, p_diff);
    drawnow limitrate

    if maxDiff < steadyTolerance
        fprintf('Steady-state reached at time %d (Tolerance = %.3e)\n', n * dt, maxDiff);
        steadyReached = true;
    end

    u_old = u;
    v_old = v;
    p_old = p;
    totalIterations = totalIterations + iterations;

end

elapsedTime = toc;  % Stop timer and get elapsed seconds
fprintf('Total elapsed time: %.2f seconds\n', elapsedTime);

%% Plot setup
% x grid
x = dx/2 : dx : length - dx/2;

% y grid (inverted)
y = height - (dy/2 : dy : height - dy/2);

[xm, ym] = meshgrid(x, y);

% Indices of interior cell centers
iC = 2:numCellsY+1;
jC = 2:numCellsX+1;

% Interpolate staggered u, v to cell centers
u_center = 0.5 * (u_new(iC, jC) + u_new(iC, jC+1));
v_center = 0.5 * (v_new(iC, jC) + v_new(iC+1, jC));
p_center = p_new(iC, jC);

%% Contour plots
% Velocity magnitude
U_center_mag = sqrt(u_center.^2 + v_center.^2);

figure(1);
contourf(x, y, U_center_mag, 20, 'LineColor', 'none'); % 20 contour levels
colorbar;
colormap(jet(21));
hold on;
h = streamslice(xm, ym, u_center, v_center, 2, 'cubic');
set(h,'Color','w');  % make streamlines white
axis square;
axis([0 1 0 1]);
xlabel('x'); ylabel('y');
title('Streamlines + Velocity Magnitude');
grid on;

% u velocity
figure(2);
contourf(x, y, u_center, 21, 'LineColor', 'none');
colorbar;
colormap(jet(21));
axis square;
xlabel('x'); ylabel('y');
title('u velocity Contour');

% v velocity
figure(3);
contourf(x, y, v_center, 21, 'LineColor', 'none');
colorbar;
colormap(jet(21));
axis square;
xlabel('x'); ylabel('y');
title('v velocity Contour');

% pressure
figure(4);
contourf(x, y, p_center, 21, 'LineColor', 'none');
colorbar;
colormap(jet(21));
axis square;
xlabel('x'); ylabel('y');
title('Pressure Contour');

%% Extract centerline velocities for comparison with Ghia et al.
% Middle indices
j_mid = round(numCellsX/2);
i_mid = round(numCellsY/2);

% Extract u vertical centerline (x = 0.5)
u_centerline = u_center(:, j_mid);
y_center = y(:);

% Extract v horizontal centerline (y = 0.5)
v_centerline = v_center(i_mid, :);
x_center = x(:);

%% Compare with Ghia
% Read Ghia et al. (1982) benchmark values 
[scriptDir, fileName, ~] = fileparts(mfilename('fullpath')); % Get file path and file name
ghiaFolder = fullfile(scriptDir, 'Ghia Centerline');
u_table = readtable(fullfile(ghiaFolder, 'Ghia-u-centerline.csv'), ReadVariableNames=true, VariableNamingRule="preserve");
v_table = readtable(fullfile(ghiaFolder, 'Ghia-v-centerline.csv'), ReadVariableNames=true, VariableNamingRule="preserve");

% Extract data
ReString = num2str(Re, '%g');
ghia_y = u_table.y;
ghia_x = v_table.x;
ghia_u = u_table.(ReString);
ghia_v = v_table.(ReString);

figure(5);
plot(u_centerline, y_center, 'b'); hold on;
plot(ghia_u, ghia_y, 'o'); % Re=100
axis square;
xlabel('u');
ylabel('y');
legend('SIMPLE', 'Ghia et al.');
title(['Vertical Centerline u Velocity Re = ', ReString]);
grid on;

figure(6);
plot(x_center, v_centerline, 'b'); hold on;
plot(ghia_x, ghia_v, 'o');
axis square;
xlabel('x');
ylabel('v');
legend('SIMPLE','Ghia et al.');
title(['Horizontal Centerline v Velocity Re = ', ReString]);
grid on;

%% Save results
% Create results folder
[~, branchName] = system('git rev-parse --abbrev-ref HEAD');
branchName = strtrim(branchName); % Get branch name

if strcmp(branchName, "master")
    branchName = "Euler";
end

filePrefix = "Re" + ReString + "_";
timeStamp = string(datetime("now", "Format", "dd-MMM_HH-mm")); % Get current time (e.g. 19Jan_02-38)

resultsFolder = fullfile(scriptDir, "Results_" + fileName + "_" + branchName, filePrefix + timeStamp);

mkdir(resultsFolder);   % Create the folder
disp(["Folder created: " resultsFolder]);

% Performance Log
logFile = fullfile(resultsFolder, filePrefix + "Performance_Log.txt");
fileID = fopen(logFile, 'w');
fprintf(fileID, 'SIMULATION PERFORMANCE LOG\n');
fprintf(fileID, '==========================\n');
fprintf(fileID, "File: %s\n", fileName + "_" + branchName);
fprintf(fileID, 'Reynolds Number: %d\n', Re);
fprintf(fileID, 'Grid Size: %d x %d\n', numCellsX, numCellsY);
fprintf(fileID, 'Courant Number: %g\n', CFL);
fprintf(fileID, 'Time Step (dt): %.1e\n', dt);
fprintf(fileID, 'Pressure Under-relaxation: %g\n', alpha_p);
fprintf(fileID, 'Total Time Steps Completed: %d\n', n);
fprintf(fileID, 'Total Time Completed: %.4f seconds\n', n * dt);
fprintf(fileID, 'Total SIMPLE Iterations: %d\n', totalIterations);
fprintf(fileID, 'Total Elapsed Time: %.2f seconds\n', elapsedTime);
fprintf(fileID, 'Final Steady State Tolerance %g\n', maxDiff);

if n > 0
    fprintf(fileID, 'Avg Iterations per Time Step: %.2f\n', totalIterations / n);
    fprintf(fileID, 'Avg Time per Time Step: %.4f seconds\n', elapsedTime / n);
end

fclose(fileID);

% Exporting results
set(groot, 'defaultAxesToolbarVisible', 'off'); % Set the default figure property for the Axes Toolbar to 'never'
writematrix(u_centerline, fullfile(resultsFolder, filePrefix + "u_Center_Line.csv"));
writematrix(v_centerline, fullfile(resultsFolder, filePrefix + "v_Center_Line.csv"));
exportgraphics(figure(1), fullfile(resultsFolder, filePrefix + "Velocity_Magnitude.png"), 'Resolution', 300);
exportgraphics(figure(2), fullfile(resultsFolder, filePrefix + "u_Velocity.png"), 'Resolution', 300);
exportgraphics(figure(3), fullfile(resultsFolder, filePrefix + "v_Velocity.png"), 'Resolution', 300);
exportgraphics(figure(4), fullfile(resultsFolder, filePrefix + "Pressure.png"), 'Resolution', 300);
exportgraphics(figure(5), fullfile(resultsFolder, filePrefix + "u_Centerline.png"), 'Resolution', 300);
exportgraphics(figure(6), fullfile(resultsFolder, filePrefix + "v_Centerline.png"), 'Resolution', 300);