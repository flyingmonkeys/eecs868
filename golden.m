function alpha_final = golden(x,d,c1,c2,alpha_max)
% Return the alpha that minimizes f() using the Golden Section line
% search method and meets strong Wolfe conditions as specified by c1 and c2
% x = starting point vector
% d = direction vector
% c1 = Armijo parameter
% c2 = Curvature parameter
% alpha_max = upper limit for alpha
% Note: Requires functions to return f(), f'(), phi'()

N = 20; % just in case... limit the max number of iterations

% Simulation parameters
a_min = 0.0;        % Given in problem
c = zeros(2,1);     % Endpoints of the Golden Section search
e = zeros(N,1);     % candidate points
alpha = zeros(N,1); % history of alpha(k) for plotting
phi = zeros(N,1);   % history of objective function phi(k) for plotting
phip = zeros(N,1);  % history of phi'(k) for plotting

% Golden section search
tau2 = 2.0/(1.0+sqrt(5.0));
c(1) = a_min;       % c1 < c2 always
c(2) = alpha_max;

% Initialize first two points
du = c(2) - c(1);       % uncertainty distance
e(1) = c(1) + tau2*du;
e(2) = c(2) - tau2*du;
lowest = e(1);          % arbitrarily assign lowest to first candidate point
y = f(x + lowest*d);
alpha(1) = c(1) + du/2;    % as defined in the project description
phi(1) = f(x + alpha(1)*d);
phip(1) = phiprime(alpha(1),d,x);

% Iterate until Wolfe condition is met
for k=2:N
	% Move in boundaries
	if( f(x + e(k)*d) < y )
		new = farthest(e(k),c);
		c(new) = lowest;        % move in boundary
		lowest = e(k);        % mark index of lowest value
		y = f(x + lowest*d);	% record the lowest value
	else
		new = farthest(lowest,c);
		c(new) = e(k);	% move in boundary
	end

	% Get next evaluation point
	du = c(2) - c(1);	% Calculate new uncertainty width
	if( new == 1 )
		e(k+1) = c(1) + tau2*du;
	else
		e(k+1) = c(2) - tau2*du;
    	end
    
    alpha(k) = c(1) + du/2;    % as defined in the project description
    phi(k) = f(x + alpha(k)*d);
    phip(k) = phiprime(alpha(k),d,x);

    % Check Wolfe conditions
    fprintf('k=%2i alpha=%6.3f f(alpha)=%6.3f, d=%6.3f: ',k,alpha(k),f(x + alpha(k)*d),du);
    % Armijo condition check
    if( f(x + alpha(k)*d) <= f(x) + c1*alpha(k)*phiprime(0,d,x) )
        fprintf('(Armijo checks)');
    end
    % Curvature condition check
 %   if( abs(phiprime(alpha(k),d,x)) <= c2*abs(phiprime(0,d,x)) )
    if( phiprime(alpha(k),d,x) >= c2*phiprime(0,d,x) )
        fprintf('(Curvature checks)');
    end
    fprintf('\n');
end

alpha_final = alpha(k);

% Plots
k = [1:1:N];
figure(1)
plot(k,alpha);
title('Alpha_k vs. Iteration');
xlabel('Iteration');
ylabel('Alpha');
figure(2)
plot(k,phi);
title('Phi_k (Alpha) vs. Iteration');
xlabel('Iteration');
ylabel('Phi(Alpha)');
figure(3)
plot(k,phip);
title('Phi''_k (Alpha) vs. Iteration');
xlabel('Iteration');
ylabel('Phi''(Alpha)');

end


