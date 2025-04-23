% Load the result for Re = 100 and N = 128
load('Re1000_N64.mat');  % must contain uf, vf, dx, dy

% Construct the non-uniform physical grid
N = length(dx) - 2;
x = zeros(1, N+2);
y = zeros(1, N+2);

for i = 2:N+2
    x(i) = x(i-1) + 0.5 * (dx(i-1) + dx(i));
    y(i) = y(i-1) + 0.5 * (dy(i-1) + dy(i));
end

[X, Y] = meshgrid(x, y);

% Define uniform grid for smooth plotting
num_points = 300;  % Higher resolution for fine grid
x_uni = linspace(0, 1, num_points);
y_uni = linspace(0, 1, num_points);
[Xq, Yq] = meshgrid(x_uni, y_uni);

% Interpolate velocity data to the uniform grid
U_interp = griddata(X, Y, uf, Xq, Yq, 'cubic');
V_interp = griddata(X, Y, vf, Xq, Yq, 'cubic');

% Plot the streamlines
figure;
streamslice(Xq, Yq, U_interp, V_interp);
axis equal;
set(gca, 'YDir', 'reverse');  % Flip Y-axis for traditional cavity view
title('Streamlines for Lid Driven Cavity (Re=1000, N=64)');
xlabel('X'); ylabel('Y');
