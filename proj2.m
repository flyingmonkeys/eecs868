% EECS 868 Project 2
% Top-level script to set up simulation parameters for executing a Steepest
% Descent algorithm using a Golden Section line search with Strong Wolfe
% stopping conditions on the Rosenbrock function.
% Line is defined within the functions f.m and phiprime.m

clear all;
close all;

pout = 0;

% Simulation parameters
c1    = 0.40;            % Wolfe parameter (Armijo)
c2    = 0.70;            % Wolfe parameter (curvature)
a_max = 1.9;            % maximum alpha value to search for
g     = zeros(2,1);     % gradient vector
%x     = [1.2 1.2]';     % initial point 0.01,0.55,2 6 iter
%x     = [-1.2 1.0]';    % initial point 0.01,0.55,2 11 iter (c2=0.35 144 iter)
%x     = [10 0]';        % initial point 0.01,0.55,2 12 iter
%x     = [1.5 15]';      % initial point 0.01,0.55,8 169 iter 0.4,0.7,2.5 12899 iter
x = [3.74125 14]';
N     = 40000;          % max number of iterations (just in case)
num_iter = 0;           % total number of low-level iterations

x_hist = zeros(N,2);    % history of x per iteration
g_hist = zeros(N,1);    % history of gradient vector length

% Get initial conditions
% Find direction vector via negative gradient
g(1) = 2 * (x(1) - 1 - 200*x(1)*(x(2)-(x(1)^2)));
g(2) = 200 * (x(2) - (x(1)^2));
d    = -g;	% direction vector

if( pout )
    fprintf('x = [%6.3f %6.3f]\n',x(1),x(2));
    fprintf('g=[%6.3f %6.3f] ||g|| = %6.3f d=[%6.3f %6.3f]\n', g(1),g(2),sqrt(g'*g),d(1)/norm(d),d(2)/norm(d));
end

for k=1:N
    % Record histories
    x_hist(k,:) = x';
    g_hist(k,:) = (g'*g);

    % Golden Section search call
    [alpha,iters] = golden(x,d/norm(d),c1,c2,a_max);
    num_iter = num_iter + iters; % track number of low-level iterations

    % Update solution
    x = x + alpha*d/norm(d);

    % Find gradient vector
    g_k  = g; % record old gradient value
    g(1) = 2 * (x(1) - 1 - 200*x(1)*(x(2)-(x(1)^2)));
    g(2) = 200 * (x(2) - (x(1)^2));

    % Check for termination condition
    if( sqrt(g'*g) <= 0.2 )
        break;
    end

    % Calculate beta
    beta = (g'*g) / (g_k'*g_k);

    % Update direction vector
    d = -g;             % steepest descent update
%	d = -g + beta*d;     % Fletcher-Reeves update

    if( pout )
        fprintf('x = [%8.5f %8.5f]\n',x(1),x(2));
        fprintf('g = [%6.3f %6.3f] ||g|| = %6.3f d=[%6.3f %6.3f]\n', g(1),g(2),sqrt(g'*g),d(1)/norm(d),d(2)/norm(d));
    end
    
end

% Final results
x_hist(k+1,:) = x';
fprintf('x    = [%8.5f %8.5f]\n',x(1),x(2));
fprintf('f(x) = %8.6f\n',f(x));
fprintf('Number of high-level iterations = %i\n',k);
fprintf('Number of low-level iterations  = %i\n',num_iter);

%% Plots
figure(1)
plot(x_hist(1:k+1,1),x_hist(1:k+1,2),'-o');
grid on;



