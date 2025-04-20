% Parameters
N = 6;                  % Number of nodes
R = 1;                  % Initial radius
F0 = 1;                 % Expansion force
gamma = 0.01;           % Surface tension coefficient
alpha = 0.5;            % Friction coefficient
c = 0.1 * ones(N, 1);   % E. coli concentration at each node
dt = 0.01;              % Time step
T = 100;                % Number of time steps

% Initial positions (circle)
theta = linspace(0, 2*pi, N+1)';
theta = theta(1:end-1); % Remove duplicate at 2*pi
x = R * cos(theta);
y = R * sin(theta);

% Simulation loop
for t = 1:T
    % Compute edges and curvatures
    dx = circshift(x, -1) - x;
    dy = circshift(y, -1) - y;
    lengths = sqrt(dx.^2 + dy.^2);
    normals_x = -dy ./ lengths;
    normals_y = dx ./ lengths;
    curvatures = abs((circshift(dx, -1) ./ circshift(lengths, -1)) - (dx ./ lengths)) ./ ...
                 (0.5 * (lengths + circshift(lengths, -1)));
    
    % Compute forces
    fx = normals_x .* (F0 - gamma * curvatures) ./ (1 + alpha * c);
    fy = normals_y .* (F0 - gamma * curvatures) ./ (1 + alpha * c);
    
    % Update positions
    x = x + dt * fx;
    y = y + dt * fy;
    
    % Visualization
    clf;
    plot([x; x(1)], [y; y(1)], '-o');
    axis equal;
    title(['Time step: ', num2str(t)]);
    drawnow;
    print(x)
end
