% Load results
load('Re100_N128.mat');     % contains: uf, vf, dx, dy
uf_100 = uf; vf_100 = vf; dx_100 = dx; dy_100 = dy;

load('Re1000_N128.mat');    % contains: uf, vf, dx, dy
uf_1000 = uf; vf_1000 = vf; dx_1000 = dx; dy_1000 = dy;

% Function to compute and plot vorticity
function plotVorticity(uf, vf, dx, dy, Re)
    N = length(dx) - 2;
    % Generate physical non-uniform grid
    x = zeros(1, N+2);
    y = zeros(1, N+2);
    for i = 2:N+2
        x(i) = x(i-1) + 0.5 * (dx(i-1) + dx(i));
        y(i) = y(i-1) + 0.5 * (dy(i-1) + dy(i));
    end
    [X, Y] = meshgrid(x, y);

    % Compute vorticity on coarse grid
    [rows, cols] = size(uf);
    omega = zeros(rows, cols);
    for i = 2:rows-1
        for j = 2:cols-1
            dvdx = (vf(i, j+1) - vf(i, j-1)) / (x(j+1) - x(j-1));
            dudy = (uf(i+1, j) - uf(i-1, j)) / (y(i+1) - y(i-1));
            omega(i, j) = dvdx - dudy;
        end
    end

    % Interpolate to a fine uniform grid for smooth plotting
    num_points = 300;
    x_uni = linspace(0, 1, num_points);
    y_uni = linspace(0, 1, num_points);
    [Xq, Yq] = meshgrid(x_uni, y_uni);
    omega_interp = griddata(X, Y, omega, Xq, Yq, 'cubic');

    % Plot vorticity contours
    figure;
    contourf(Xq, Yq, omega_interp, 30, 'LineColor', 'none');
    colorbar;
    axis equal tight;
    set(gca, 'YDir', 'reverse');
    xlabel('X'); ylabel('Y');
    title(['Vorticity Contours, Re = ', num2str(Re), ', N = ', num2str(N)]);
end

% Plot for Re = 100
plotVorticity(uf_100, vf_100, dx_100, dy_100, 100);

% Plot for Re = 1000
plotVorticity(uf_1000, vf_1000, dx_1000, dy_1000, 1000);
