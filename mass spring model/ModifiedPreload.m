
xp = [0, 0, 0]';
% xn0 = -10;
gxp = [0, 0, 0]';
params.func.static.preload = struct('xe0', xe0, 'Rx', Rx, 'xp', xp, 'gxp', gxp);
params.func.fc.mu = 0.3 * [1; 1]; 