function alpha_final = golden(x,d,c1,c2,alpha_max)
% Return the alpha that minimizes f() using the Golden Section line
% search method and meets strong Wolfe conditions as specified by c1 and c2
% x = starting point vector
% d = direction vector
% c1 = Armijo parameter
% c2 = Curvature parameter
% alpha_max = upper limit for alpha
% Note: Requires functions to return f(x) and phi'(x)

N = 20; % just in case... limit the max number of iterations

% Simulation parameters
a_min = 0.0;        % Given in problem
c = zeros(2,1);     % Endpoints of the Golden Section search
e = zeros(N,1);     % candidate points
alpha = zeros(N,1); % history of alpha(k) for plotting
phi = zeros(N,1);   % history of objective function phi(k) for plotting
phip = zeros(N,1);  % history of phi'(k) for plotting

% Golden section search
tau2 = 2.0/(1.0+sqrt(5.0)); % = 0.618
c(1) = a_min;               % c1 < c2 always
c(2) = alpha_max;
fprintf('c1 = %4.2f c2 = %4.2f\n',c1,c2);

% Initialize first two points
k        = 1;               % iteration number
du       = c(2) - c(1);     % uncertainty distance
e(1)     = c(1) + tau2*du;
e(2)     = c(2) - tau2*du;
lowest   = e(1);            % arbitrarily assign lowest to first candidate point
y        = f(x + lowest*d);
alpha(1) = c(1) + du/2;     % as defined in the project description
phi(1)   = f(x + alpha(1)*d);
phip(1)  = phiprime(alpha(1),d,x);
fprintf('k=%2i alpha=%6.3f f(alpha)=%6.3f, d=%6.3f: \n',k,alpha(k),f(x + alpha(k)*d),du);

% Iterate until Wolfe condition is met
for k=2:N
	% Move in boundaries
	if( f(x + e(k)*d) < y )
		new = farthest(e(k),c);
		c(new) = lowest;        % move in boundary
		lowest = e(k);          % mark index of lowest value
		y = f(x + lowest*d);	% record the lowest value
	else
		new = farthest(lowest,c);
		c(new) = e(k);          % move in boundary
	end

	% Get next evaluation point
	du = c(2) - c(1);           % Calculate new uncertainty width
    if( new == 1 )
        e(k+1) = c(1) + tau2*du;
    else
        e(k+1) = c(2) - tau2*du;
    end
    
    alpha(k) = c(1) + du/2;     % as defined in the project description
    phi(k)   = f(x + alpha(k)*d);
    phip(k)  = phiprime(alpha(k),d,x);

    % Check Wolfe conditions
    fprintf('k=%2i alpha=%6.3f f(alpha)=%6.3f, d=%6.3f: ',k,alpha(k),f(x + alpha(k)*d),du);
    checkcount = 0;
    % Armijo condition check
    if( f(x + alpha(k)*d) <= f(x) + c1*alpha(k)*phiprime(0,d,x) )
        fprintf('(Armijo checks)');
        checkcount = checkcount + 1;
    end
    % Curvature condition check
    if( abs(phiprime(alpha(k),d,x)) <= c2*abs(phiprime(0,d,x)) ) % Strong Wolfe
%    if( phiprime(alpha(k),d,x) >= c2*phiprime(0,d,x) ) % Normal Wolfe
        fprintf('(Curvature checks)');
        checkcount = checkcount + 1;
    end
    fprintf('\n');

    if( checkcount == 2 )
        fprintf('Wolfe conditions have been met!\n');
        break;
    end
end

if( k >= N )
    fprintf('Wolfe conditions not satisfied!\n');
end

alpha_final = alpha(k);

% Plots
iter = (1:1:k);

% Alpha(k) vs. k
figure(1)
plot(iter,alpha(1:k));
title('\alpha_k vs. Iteration');
xlabel('Iteration');
ylabel('\alpha');
grid on;

% Phi(alpha(k)) vs. k
figure(2)
plot(iter,phi(1:k));
title('\phi_k (\alpha) vs. Iteration');
xlabel('Iteration');
ylabel('\phi(\alpha)');
grid on;

% Phi'(alpha(k)) vs. k
figure(3)
plot(iter,phip(1:k));
title('\phi''_k (\alpha) vs. Iteration');
xlabel('Iteration');
ylabel('\phi''(\alpha)');
grid on;

end


