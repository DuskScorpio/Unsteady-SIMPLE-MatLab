function PlotStaggeredGrid()
% PlotStaggeredGrid  Visualises the staggered-grid layout for a 3x3
% lid-driven cavity simulation, including ghost cells.
%
% Array sizes match LidDrivenCavity.m with numCellsX = numCellsY = 3:
%   u  : (Ny+2) x (Nx+3)   faces on vertical (east/west) cell boundaries
%   v  : (Ny+3) x (Nx+2)   faces on horizontal (north/south) cell boundaries
%   p  : (Ny+2) x (Nx+2)   cell centres (including 1 ghost layer)

    %% ── Grid parameters ──────────────────────────────────────────────────
    Nx = 3;          % interior cells in x
    Ny = 3;          % interior cells in y
    dx = 1;          % cell width  (unit grid for clarity)
    dy = 1;          % cell height

    % Physical domain corners (interior cells only)
    x0 = 0;          % left   boundary
    y0 = 0;          % bottom boundary
    xL = Nx * dx;    % right  boundary
    yL = Ny * dy;    % top    boundary

    % ── Cell-centre coordinates ───────────────────────────────────────────
    % Interior cell centres  (indices 2..Ny+1 in the arrays)
    xC_int = x0 + (0.5:Nx) * dx;          % 1 x Nx
    yC_int = yL - (0.5:Ny) * dy;          % 1 x Ny  (row 2 = top interior row)

    % Ghost cell centres (one layer each side)
    xC_ghost_W = x0 - 0.5*dx;
    xC_ghost_E = xL + 0.5*dx;
    yC_ghost_N = yL + 0.5*dy;
    yC_ghost_S = y0 - 0.5*dy;

    % All cell centres in row/column order (ghost + interior + ghost)
    xC_all = [xC_ghost_W, xC_int, xC_ghost_E];   % 1 x (Nx+2)
    yC_all = [yC_ghost_N, yC_int, yC_ghost_S];   % 1 x (Ny+2)  top→bottom

    % ── u-node coordinates ────────────────────────────────────────────────
    % u lives on east/west cell faces.
    % Array columns run from j=1..Nx+3:
    %   j=1        : extra ghost (west-west face)
    %   j=2        : west ghost face  (x = x0)
    %   j=3..Nx+2  : interior faces   (x = x0+dx .. x0+Nx*dx-dx, i.e. between cells)
    %                Actually interior u is at x0+dx, x0+2dx ... for Nx=3 → x=1,2
    %                plus the right boundary face at x0+Nx*dx = 3
    %   j=Nx+3     : east ghost face
    %
    % Easier: u column j sits at x = x0 + (j-2)*dx   (j=1..Nx+3)
    xu_cols = x0 + ((1:Nx+3) - 2) * dx;      % 1 x (Nx+3)
    % u row i sits at the same y as pressure row i
    yu_rows = yC_all;                          % 1 x (Ny+2)

    % ── v-node coordinates ────────────────────────────────────────────────
    % v lives on north/south cell faces.
    % Array rows run from i=1..Ny+3:
    %   i=1       : extra ghost (north-north face)
    %   i=2       : north ghost face (y = yL)
    %   i=3..Ny+2 : interior faces
    %   i=Ny+3    : south ghost face
    %
    % v row i sits at y = yL - (i-2)*dy
    yv_rows = yL - ((1:Ny+3) - 2) * dy;       % 1 x (Ny+3)
    xv_cols = xC_all;                          % 1 x (Nx+2)

    %% ── Figure setup ─────────────────────────────────────────────────────
    fig = figure('Name','Staggered Grid – 3x3 with Ghost Cells', ...
                 'Units','normalized','OuterPosition',[0.05 0.05 0.9 0.9]);
    fig.Theme = 'light';
    hold on; axis equal;

    %% ── Draw cell boundaries (dotted lines only, interior cells) ─────────
    % Clearly visible gray lines with good contrast
    lineOpts = {'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',1.2};

    % Vertical lines spanning interior domain height
    for j = 0:Nx
        xx = x0 + j*dx;
        plot([xx xx],[y0 yL], lineOpts{:});
    end
    % Horizontal lines spanning interior domain width
    for i = 0:Ny
        yy = y0 + i*dy;
        plot([x0 xL],[yy yy], lineOpts{:});
    end

    %% ── Shade interior cells only ────────────────────────────────────────
    for ci = 1:Nx
        for cj = 1:Ny
            xPatch = x0 + (ci-1)*dx;
            yPatch = y0 + (cj-1)*dy;
            patch([xPatch xPatch+dx xPatch+dx xPatch], ...
                  [yPatch yPatch yPatch+dy yPatch+dy], ...
                  [0.94 0.96 1.0], 'FaceAlpha', 0.55, 'EdgeColor','none');
        end
    end
    %% ── Helper: arrow ────────────────────────────────────────────────────
    arrowLen_u = 0.35 * dx;   % arrow half-length for u (horizontal)
    arrowLen_v = 0.35 * dy;   % arrow half-length for v (vertical)

    %% ── Plot u nodes (blue horizontal arrows) ────────────────────────────
    % u array is (Ny+2) x (Nx+3).  Row i ↔ yC_all(i), col j ↔ xu_cols(j).
    for i = 1:Ny+2
        for j = 1:Nx+3
            x_u = xu_cols(j);
            y_u = yu_rows(i);
            isGhost = (i==1)||(i==Ny+2)||(j==1)||(j==Nx+3);
            col = [0.00 0.45 1.00];     % vivid blue
            if isGhost, col = col * 0.4 + 0.6; end  % lighter for ghost
            quiver(x_u - arrowLen_u, y_u, 2*arrowLen_u, 0, 0, ...
                   'Color',col,'LineWidth',1.0,'MaxHeadSize',0.6);
            label = sprintf('u(%d,%d)',i,j);
            text(x_u, y_u + 0.18*dy, label, ...
                 'FontSize',8,'HorizontalAlignment','center', ...
                 'Color',col,'FontWeight','normal');
        end
    end

    %% ── Plot v nodes (green vertical arrows) ─────────────────────────────
    % v array is (Ny+3) x (Nx+2).  Row i ↔ yv_rows(i), col j ↔ xv_cols(j).
    for i = 1:Ny+3
        for j = 1:Nx+2
            x_v = xv_cols(j);
            y_v = yv_rows(i);
            isGhost = (i==1)||(i==Ny+3)||(j==1)||(j==Nx+2);
            col = [0.00 0.72 0.10];     % vivid green
            if isGhost, col = col * 0.35 + 0.65; end
            quiver(x_v, y_v + arrowLen_v, 0, -2*arrowLen_v, 0, ...
                   'Color',col,'LineWidth',1.0,'MaxHeadSize',0.6);
            label = sprintf('v(%d,%d)',i,j);
            text(x_v - 0.20*dx, y_v, label, ...
                 'FontSize',8,'HorizontalAlignment','right', ...
                 'Color',col,'FontWeight','normal');
        end
    end

    %% ── Plot p nodes (red dots) ──────────────────────────────────────────
    % p array is (Ny+2) x (Nx+2).  Row i ↔ yC_all(i), col j ↔ xC_all(j).
    for i = 1:Ny+2
        for j = 1:Nx+2
            x_p = xC_all(j);
            y_p = yC_all(i);
            isGhost = (i==1)||(i==Ny+2)||(j==1)||(j==Nx+2);
            mSize = 8;
            if isGhost
                plot(x_p, y_p, 'o', 'MarkerSize', mSize, ...
                     'MarkerFaceColor',[1 0.6 0.6], ...
                     'MarkerEdgeColor',[0.8 0.2 0.2],'LineWidth',0.8);
            else
                plot(x_p, y_p, 'o', 'MarkerSize', mSize, ...
                     'MarkerFaceColor',[0.85 0 0], ...
                     'MarkerEdgeColor',[0.5 0 0],'LineWidth',1.0);
            end
            label = sprintf('P(%d,%d)',i,j);
            text(x_p + 0.04*dx, y_p - 0.18*dy, label, ...
                 'FontSize',8,'HorizontalAlignment','left', ...
                 'Color',[0.7 0 0]);
        end
    end

    %% ── Add i,j axis on top-left corner using annotations ────────────────
    % Position axis OUTSIDE the plot frame with equal horizontal and vertical gaps
    % j-axis arrow (horizontal, blue)
    annotation('arrow', [0.305 0.348], [0.960 0.960], ...
               'Color', [0 0.45 1.00], 'LineWidth', 1.0, 'HeadLength', 4, 'HeadWidth', 4);
    annotation('textbox', [0.354 0.953 0.02 0.02], 'String', 'j', ...
               'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
               'Color', [0 0.45 1.00], 'HorizontalAlignment', 'left');
    
    % i-axis arrow (vertical downward, green) - adjusted for visual equality
    annotation('arrow', [0.305 0.305], [0.960 0.902], ...
               'Color', [0.00 0.72 0.10], 'LineWidth', 1.0, 'HeadLength', 4, 'HeadWidth', 4);
    annotation('textbox', [0.298 0.895 0.02 0.02], 'String', 'i', ...
               'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
               'Color', [0.00 0.72 0.10], 'HorizontalAlignment', 'center');
    
    %% ── Axes, labels, legend ─────────────────────────────────────────────
    xLim = [x0 - 1.6*dx,  xL + 1.6*dx];
    yLim = [y0 - 1.6*dy,  yL + 1.6*dy];
    xlim(xLim);  ylim(yLim);
    xlabel('x (cell index)','FontSize',14);
    ylabel('y (cell index)','FontSize',14);

    % Tick labels in cell-index units
    xticks(x0-dx : dx : xL+dx);
    yticks(y0-dy : dy : yL+dy);

    % Legend proxies with clear colored lines
    h1 = plot(nan,nan,'--','Color',[0.4 0.4 0.4],'LineWidth',1.5);
    h2 = plot(nan,nan,'o','MarkerFaceColor',[0.85 0 0], ...
              'MarkerEdgeColor',[0.5 0 0],'MarkerSize',7,'LineWidth',1);
    h3 = plot(nan,nan,'-','Color',[0.00 0.45 1.00],'LineWidth',2.5);
    h4 = plot(nan,nan,'-','Color',[0.00 0.72 0.10],'LineWidth',2.5);
    lgd = legend([h1 h2 h3 h4], ...
           {'Cell boundary','Pressure (P)', ...
            'x-velocity (u)','y-velocity (v)'}, ...
           'Location','southoutside','Orientation','horizontal','FontSize',12);
    lgd.Box = 'off';

    hold off;
    grid off;
    box on;

    % Export without figure chrome and toolbar
    outFile = fullfile(fileparts(mfilename('fullpath')), 'StaggeredGrid_3x3.png');
    exportgraphics(gca, outFile, 'Resolution', 300, 'BackgroundColor', 'white', ...
                   'ContentType', 'vector');
    fprintf('Figure saved to: %s\n', outFile);
end