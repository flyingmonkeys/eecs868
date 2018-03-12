function out = f(x)
% Return the scalar value of the objective function
% f(x1,x2) = exp((x2-x1)/2) + (5/8)(x1-x2)

%proj1	out = exp((x(2)-x(1))/2) + (5/8)*(x(1)-x(2));
	out = 100*(x(2)-(x(1)^2))^2 + (1-x(1))^2;
end


