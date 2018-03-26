%close all;
hfig = figure(1);

s = 0.05;
x_min = -10;
x_max = 10;
y_min = -10;
y_max = 20;
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

%surf(X, Y, Z, C, 'EdgeColor', 'none', 'LineStyle', 'none');
v = [0 0.1 0.2 0.3 0.4 0.5 1 1.001 1.002 1.01 1.02 1.03 1.1 1.5 2.0 3.0 4.0 5.0 10 20 30 40 50 100 200 500 1000 2000 5000 10000 50000 200000];
%contour(X,Y,Z,40);
contour(X,Y,Z,v);

grid on;
title('Rosenbrock Contour');
colormap = jet;
 
%axis([x_min, x_max, y_min, y_max, 0, maxZ]);
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

%% Plots
%
hold on;
plot(x_hist(1:k+1,1),x_hist(1:k+1,2),'-o','Color','r');
grid on;
%}