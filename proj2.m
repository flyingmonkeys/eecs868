% EECS 868 Project 2
% Top-level script to set up simulation parameters for executing a Golden
% Section line search with Strong Wolfe stopping conditions.
% Line is defined within the functions f.m and phiprime.m

clear all;
close all;

% Simulation parameters
c1 = 0.4;                      % Wolfe parameter (Armijo)
c2 = 0.7;                      % Wolfe parameter (curvature)
%c2 = 0.7;
%a_max = 0.85;                   % maximum alpha value to search for
a_max = 2.5;
g = zeros(2,1);			% gradient vector
%x = [1.2 1.2]';                 % initial point 0.01,0.55,2 6 iter
%x = [-1.2 1.0]';                % initial point 0.01,0.55,2 11 iter (c2=0.35 144 iter)
%x = [10 0]';                    % initial point 0.01,0.55,2 12 iter
x = [1.5 15]';                  % initial point 0.01,0.55,8 169 iter 0.4,0.7,2.5 12899 iter
N = 30000;			% max number of iterations (just in case)

% Get initial conditions
% Find direction vector via negative gradient
g(1) = 2 * (x(1) - 1 - 200*x(1)*(x(2)-(x(1)^2)));
g(2) = 200 * (x(2) - (x(1)^2));
d    = -g;	% direction vector

%fprintf('x = [%6.3f %6.3f]\n',x(1),x(2));
%fprintf('g=[%6.3f %6.3f] ||g|| = %6.3f d=[%6.3f %6.3f]\n', g(1),g(2),sqrt(g'*g),d(1)/norm(d),d(2)/norm(d));

for k=1:N
	% Golden Section search call
	alpha = golden(x,d/norm(d),c1,c2,a_max);
%	fprintf('Final alpha at search termination = %6.3f\n',alpha);

	% Update solution
	x = x + alpha*d/norm(d);

	% Find gradient vector
	g_k  = g; % record old gradient value
	g(1) = 2 * (x(1) - 1 - 200*x(1)*(x(2)-(x(1)^2)));
	g(2) = 200 * (x(2) - (x(1)^2));

%	fprintf('x = [%8.5f %8.5f]\n',x(1),x(2));

	% Check for termination condition
	if( sqrt(g'*g) <= 0.2 )
%		g
		fprintf('Algorithm termination after %i iterations!\n',k);
		break;
	end

	% Calculate beta
	beta = (g'*g) / (g_k'*g_k);

	% Update direction vector
	d = -g;			% steepest descent update
%	d = -g + beta*d;	% Fletcher-Reeves update

%	fprintf('g=[%6.3f %6.3f] ||g|| = %6.3f d=[%6.3f %6.3f]\n', g(1),g(2),sqrt(g'*g),d(1)/norm(d),d(2)/norm(d));

end




