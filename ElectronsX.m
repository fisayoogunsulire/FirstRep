clearvars, clc

% Configuration Variables

electron_num = 5;

Mu = 1;

tspan = 0 : 0.00001 : 20;

% Initial State Calculation

theta0 = 360 * rand(electron_num,1);
phi0 = 360 * rand(electron_num,1);

rand_position = [sind(phi0) .* cosd(theta0), sind(phi0) .* sind(theta0), cosd(phi0)];

velocity = zeros([electron_num,3]);

initial_block = [rand_position, velocity];

initial_state = reshape(initial_block,electron_num * 6,1);

% Solve

[t, state] = ode45(@(t,state) electron_derivatives(t, state, electron_num, Mu), tspan, initial_state);

% Axis Limits

lims = [-2 2];

% Setup

figure('Color', 'w', 'Position', [100, 100, 1000, 800])

% Colors 

vib_colors = ["r";...
              "g";...
              "b";...
              "y";...
              "c";...
              "m";...
              "w"];
% Sphere Creation

[X,Y,Z] = sphere(20);
hold on
sphere_size = 0.1;
sphere_b = zeros(electron_num,1);
for j = 1:electron_num
    hold on
    sphere_b(j) = surf(X*sphere_size + state(1,0+j),...
                     Y*sphere_size + state(1,electron_num+j),...
                     Z*sphere_size + state(1,2*electron_num+j),...
                     'FaceColor', vib_colors(1+mod(j,7),:), 'EdgeColor', 'none', 'FaceAlpha', 0.9);
end


% Axes

xlabel('X Position', 'FontSize', 12)
ylabel('Y Position', 'FontSize', 12)
zlabel('Z Position', 'FontSize', 12)
grid on
axis equal
xlim(lims)
ylim(lims)
zlim(lims)
lighting gouraud
light('Position', [1, 1, 1])

% Animation

for i = 1:10000:length(t)
    
    for j = 1:electron_num

        set(sphere_b(j), 'XData', X*sphere_size + state(i,j), 'YData', Y*sphere_size + state(i,electron_num+j), 'ZData', Z*sphere_size + state(i,2*electron_num+j));
    end 
    % Update view angle
    view_angle = 30 + t(i) * 10;
    view(view_angle, 20)
    
    % Update title
    title(sprintf('Thompson Problem | Time: %.2f s', t(i)), 'FontSize', 14)
    
    drawnow;

    pause(0.1); % Control the animation speed
end

hold off