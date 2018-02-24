function alpha = golden(x,d,c1,c2,alpha_max)
% Return the alpha that minimizes f() using the Golden Section line
% search method and meets strong Wolfe conditions as specified by c1 and c2

N = 20; % just in case... limit the max number of iterations

% Simulation parameters
a_min = 0.0;
c = zeros(2,1);     % Endpoints of the Golden Section search
e = zeros(N,1);     % candidate points

% Golden section search
tau2 = 0.618;
c(1) = a_min;       % c1 < c2 always
c(2) = alpha_max;

% Initialize first two points
du = c(2) - c(1);       % uncertainty distance
e(1) = c(1) + tau2*du;
e(2) = c(2) - tau2*du;
lowest = e(1);          % arbitrarily assign lowest to first candidate point
y = f(x + lowest*d);

% Iterate until Wolfe condition is met
for k=3:N+1
	% Move in boundaries
	if( f(x + e(k-1)*d) < y )
		new = farthest(e(k-1),c);
		c(new) = lowest;        % move in boundary
		lowest = e(k-1);        % mark index of lowest value
		y = f(x + lowest*d);	% record the lowest value

	else
		new = farthest(lowest,c);
		c(new) = e(k-1);	% move in boundary
	end

	% Get next evaluation point
	du = c(2) - c(1);
	if( new == 1 )
		e(k) = c(1) + tau2*du;
	else
		e(k) = c(2) - tau2*du;
    end
    
    alpha = c(1) + du/2;    % as defined in the project description
        
    % Check Wolfe conditions
    fprintf('k=%2i alpha=%5.3f f(alpha)=%5.3f: ',k,alpha,f(x + alpha*d));
    % Armijo condition check
    if( f(x + alpha*d) <= f(x) + c1*alpha*fp(x)*d )
        fprintf('(Armijo checks)');
    end
    % Curvature condition check
    if( phiprime(alpha,d,x) <= c2*abs(phiprime(0,d,x)) )
        fprintf('(Curvature checks)');
    end
    fprintf('\n');
end

end


