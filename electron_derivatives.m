function dstate = electron_derivatives(t, state, electron_num, Mu)


% Extract values

pstate = state(1:3*electron_num);
vstate = state(3*electron_num+1:6*electron_num);

p = reshape(pstate, [electron_num,3,1]);
v = reshape(vstate, [electron_num,3,1]);

% Parameters
    
    
    K = 1;
    Q = ones(electron_num,1);
    M = ones(electron_num,1);
   
% Functions

function F = forceBtw(cord)
    
    if cord(1) == cord(2)
        F = [0,0,0];
    else
        F = (K*Q(cord(1))*Q(cord(2)) * ((p(cord(1),:)-p(cord(2),:)) /(norm(p(cord(1),:)-p(cord(2),:)))^3));
    end

end


% Displays

sum_p = norm(sum(p,1))

% Create Position Coords

[X, Y] = meshgrid(1:electron_num);

points_x = X(:);
points_y = Y(:);
points_matrix = [points_x, points_y];

C_points = num2cell(points_matrix, 2);
coords = reshape(C_points,[electron_num,electron_num]);

% CELL FUN!!!

f_cell = cellfun(@forceBtw,coords,'UniformOutput',false);
f_3d = cell2mat(reshape(f_cell, electron_num, 1,[]));
f_sq = sum(f_3d,1);
f = reshape(f_sq,3,electron_num)';

magp = vecnorm(p, 2, 2);
magv = vecnorm(v, 2, 2);
f_proj = ((dot(f,p,2) ./ magp.^2) .* p);
f_centri = ((M .* magv.^2 ./ magp.^2) .* p);
f_drag = Mu * v;
f_tot = f - f_proj - f_centri - f_drag;
a = f_tot ./ M;
dpdt = reshape(v,electron_num*3,1);
dvdt = reshape(a,electron_num*3,1);
dstate = [dpdt;dvdt];

end