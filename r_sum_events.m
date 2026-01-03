function [value,isterminal,direction] = r_sum_events(~,state,r_sum_limit, electron_num)

pstate = state(1:3*electron_num);

p = reshape(pstate, [electron_num,3,1]);

sum_r = norm(sum(p,1));

value = sum_r - r_sum_limit;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Any Direction
end