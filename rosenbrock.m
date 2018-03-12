close all;
hfig = figure(1);

s = 0.05;
x_min = 0;
x_max = 2;
y_min = 0;
y_max = 2;
X = [x_min : s : x_max+s];
Y = [y_min : s : y_max+s];
XC = X;
YC = Y;
[X, Y] = meshgrid(X, Y);
Z = (1-X).^2 + 100*(Y-X.^2).^2;

% Use log scale of Z for a batter usage of the color spectrum  
minZ = min(Z(:));  
maxZ = max(Z(:));
C = minZ + (maxZ-minZ).*log(1+Z-minZ)./log(1+maxZ-minZ); 
% See
% http://stackoverflow.com/questions/5073865/how-to-color-surface-with-stronger-contrast
% for how to color a surface with a even stronger contrast.   
% The method is as follows:
%C = Z;
%[~, index] = sort(C(:));
%C(index) = 1 : numel(index);

surf(X, Y, Z, C, 'EdgeColor', 'none', 'LineStyle', 'none');
%contour(X,Y,Z,40);
grid on;
colormap = jet;
 
axis([x_min, x_max, y_min, y_max, 0, maxZ]);
xlabel('x', 'fontsize', 18);
ylabel('y', 'fontsize', 18);
zlabel('f', 'fontsize', 18);
 
% Note that the `-dsvg' option is only supported for Simulink systems
print(hfig, '-dsvg', 'rosenbrock');
% To produce eps and pdf, use the following code. Notice that `epstopdf' may not work on Windows. 
% With octave 3.8.1, the pdf does not look very nice --- it is full of "scratches". 
%print(hfig, '-depsc', 'rosenbrock');
%system('epstopdf rosenbrock.eps');

%{
[fx,fy] = gradient(Z,0.05);
x0 = 10;
y0 = 0;
t = (XC == x0) & (YC == y0);
indt = find(t);
f_grad = [fx(indt) fy(indt)];
%}