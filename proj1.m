% EECS 868 Project 1
% Top-level script to set up simulation parameters for executing a Golden
% Section line search with Strong Wolfe stopping conditions.
% Line is defined within the functions f.m and phiprime.m

clear all;
close all;

% Simulation parameters
c1 = 0.01;                      % Wolfe parameter (Armijo)
c2 = 0.10;                      % Wolfe parameter (curvature)
a_max = 50.0;                   % maximum alpha value to search for
d = [1/sqrt(2) -1/sqrt(2)]';    % direction vector
x = [1 3]';                     % initial point

% Golden Section search call
alpha = golden(x,d,c1,c2,a_max);
fprintf('Final alpha at search termination = %6.3f\n',alpha);



