function out = phiprime(a,d,x)
% Return the scalar value of the derivative of the objective function at (x+alpha*d)
% phi'(alpha) = (d2-d1)/2 * exp((x2-x1)/2)exp((d2-d1)*a/2) - (5/4)

    out = ((d(2)-d(1))/2) * ( (exp( (x(2)-x(1))/2 )*exp( (d(2)-d(1))*a/2 )) - (5/4) );
    
end

