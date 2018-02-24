function out = fp(x)
% Return the row vector of the gradient of the objective function
% as specified in the f.m file.
% f(x1,x2) = exp((x2-x1)/2) + (5/8)(x1-x2)
% f'(x1,x2) = [-0.25(x2-x1)exp((x2-x1)/2)+(5/8) 0.25(x2-x1)exp((x2-x1)/2)-(5/8)]

    out(1) = -0.25 * (x(2) - x(1)) * exp((x(2)-x(1))/2) + (5/8);
    out(2) =  0.25 * (x(2) - x(1)) * exp((x(2)-x(1))/2) - (5/8);

end