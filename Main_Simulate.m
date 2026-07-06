% Parameters from the info
p = 1.15;
g = 0.66;
z = 0.24;
c0 = 210;        % seconds
tau0 = 1e4;      % seconds
b_bg = 1.08;     % background GR b-value
Mc = 1.5;        % lower magnitude threshold
Mmax = 8.0;      % upper magnitude cutoff
dM = 0.01;       % magnitude resolution
T = 1* 365 * 24 * 3600;  % Catalog duration in seconds

% Background rate in events per second
mu=1e-4;
seed=[];
Cata=SSAR_Simulate(T, mu, Mc, Mmax, b_bg, p, g, z, c0, tau0, dM, seed);


fprintf('Catalog Info: \n');
fprintf('Total events: %d\n', size(Cata, 1));
fprintf('Background events: %d\n', sum(Cata(:, 8) == -2));
fprintf('Aftershocks: %d\n', sum(Cata(:, 8) ~= -2));