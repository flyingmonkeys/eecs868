function [alpha_final,num_iter] = golden(x,d,c1,c2,alpha_max)
% Return the alpha that minimizes f() using the Golden Section line
% search method and meets strong Wolfe conditions as specified by c1 and c2
% Also returns the number of iterations used to to find the minimum
% x         = starting point vector
% d         = direction vector
% c1        = Armijo parameter
% c2        = Curvature parameter
% alpha_max = upper limit for alpha
% Note: Requires external functions to return f(x) and phi'(x)

pout = 0;

% Simulation parameters
N           = 200;        % just in case... limit the max number of iterations
alpha_min   = 0.0;        % Given in problem
c           = zeros(2,1); % Endpoints of the Golden Section search
e           = zeros(N,1); % Candidate points (alpha) of G.S search
alpha       = zeros(N,1); % history of alpha(k) for plotting
phi         = zeros(N,1); % history of objective function phi(k) for plotting
phip        = zeros(N,1); % history of phi'(k) for plotting

%% ------Golden section search
tau2 = 2.0/(1.0+sqrt(5.0)); % = 0.618
c(1) = alpha_min;           % c1 < c2 always
c(2) = alpha_max;
%fprintf('c1 = %4.2f c2 = %4.2f\n',c1,c2);   % for output analysis

% Initialize first two points
k        = 1;               % iteration number
du       = c(2) - c(1);     % uncertainty distance
e(1)     = c(1) + tau2*du;  % e1 = c1 + tau*d
e(2)     = c(2) - tau2*du;  % e2 = c2 - tau*d
lowest   = e(1);            % arbitrarily assign lowest alpha to first candidate point
y        = f(x + lowest*d); % y = phi(alpha)
alpha(1) = c(1) + du/2;     % as defined in the project description
phi(1)   = f(x + alpha(1)*d);
phip(1)  = phiprime(alpha(1),d,x);
if( pout )
    fprintf('k=%2i alpha=%6.3f f(alpha)=%6.3f, d=%6.3f: \n',k,alpha(k),f(x + alpha(k)*d),du);
end

% Iterate until Wolfe conditions are met (stop at N to limit iterations in
% the event of non-convergence)
for k=2:N
	% Move in boundaries
	if( f(x + e(k)*d) < y )     % if this new candidate point is lower than lowest point
		new = farthest(e(k),c); % next point will be the one farthest from this endpoint
		c(new) = lowest;        % move in boundary to the previous lowest value
		lowest = e(k);          % mark new index (alpha) of the alpha with the lowest value
		y = f(x + lowest*d);	% record the lowest value
	else
		new = farthest(lowest,c); % next point will be the one farthest from this endpoint
		c(new) = e(k);          % move in boundary
	end

	% Get next evaluation point
	du = c(2) - c(1);           % Calculate new uncertainty width
    if( new == 1 )
        e(k+1) = c(1) + tau2*du; % Next candidate point is by the lower side of the search line
    else
        e(k+1) = c(2) - tau2*du; % Next candidate point is by the higher side of the search line
    end
    
    % Log alpha, phi, and phi' for later analysis and plotting
    alpha(k) = c(1) + du/2;     % as defined in the project description
    phi(k)   = f(x + alpha(k)*d);
    phip(k)  = phiprime(alpha(k),d,x);

    % Check Wolfe stopping conditions
    if( pout )
        fprintf('k=%2i alpha=%6.3f f(alpha)=%6.3f, d=%6.3f: ',k,alpha(k),f(x + alpha(k)*d),du);
    end
    checkcount = 0; % count how many constraints have been met (0, 1, or 2)
    % Armijo condition check
    if( f(x + alpha(k)*d) <= f(x) + c1*alpha(k)*phiprime(0,d,x) )
        if( pout )
            fprintf('(Armijo checks)');
        end
        checkcount = checkcount + 1; % record that one constraint was met
    end
    % Curvature condition check (comment out either Strong Wolfe or Normal
    % Wolfe condition to switch between constraints)
    if( abs(phiprime(alpha(k),d,x)) <= c2*abs(phiprime(0,d,x)) ) % Strong Wolfe
%    if( phiprime(alpha(k),d,x) >= c2*phiprime(0,d,x) ) % Normal Wolfe
        if( pout )
            fprintf('(Curvature checks)');
        end
        checkcount = checkcount + 1; % record that one constraint was met
    end
    if( pout )
        fprintf('\n');
    end

    if( checkcount == 2 )
        if( pout )
            fprintf('Wolfe conditions have been met!\n');
        end
        break;
    end
end

if( k >= N )
    fprintf('Wolfe conditions not satisfied!!!!!!!!!!!!!!!!\n');
    fprintf('x = [%8.5f %8.5f] ',x(1),x(2));
    fprintf('d = [%6.3f %6.3f] ',d(1),d(2));
    fprintf('a=%8.6f\n',alpha(k));
end

% Mark the final value of alpha
alpha_final = alpha(k);

%% ------------------------------Plots
%{
iter = (0:1:k);     % x axis label

% Alpha(k) vs. k
figure(1)
plot(iter,[alpha_min alpha(1:k)']);
title('\alpha_k vs. Iteration');
xlabel('Iteration');
ylabel('\alpha');
grid on;

% Phi(alpha(k)) vs. k
figure(2)
plot(iter,[f(x + alpha_min*d) phi(1:k)']);
title('\phi_k (\alpha) vs. Iteration');
xlabel('Iteration');
ylabel('\phi(\alpha)');
grid on;

% Phi'(alpha(k)) vs. k
figure(3)
plot(iter,[phiprime(alpha_min,d,x) phip(1:k)']);
title('\phi''_k (\alpha) vs. Iteration');
xlabel('Iteration');
ylabel('\phi''(\alpha)');
grid on;
%}

if( nargout > 1 )
    num_iter = k;
end

end


