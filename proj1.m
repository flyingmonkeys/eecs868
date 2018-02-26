clear all;
close all;

% Simulation parameters
c1 = 0.6;                      % Wolfe parameter (Armijo)
c2 = 0.7;                       % Wolfe parameter (curvature)
a_max = 50.0;
d = [1/sqrt(2) -1/sqrt(2)]';    % direction vector
x = [1 3]';                     % initial point

% Golden search
alpha = golden(x,d,c1,c2,a_max)



