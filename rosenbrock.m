close all;
hfig = figure(1);

s = 0.05;
X = [-2 : s : 2+s];
Y = [-1 : s : 3+s];
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
colormap = jet;
 
axis([-2, 2, -1, 3, 0, 2500]);
xlabel('x', 'fontsize', 18);
ylabel('y', 'fontsize', 18);
zlabel('f', 'fontsize', 18);
 
% Note that the `-dsvg' option is only supported for Simulink systems
print(hfig, '-dsvg', 'rosenbrock');
% To produce eps and pdf, use the following code. Notice that `epstopdf' may not work on Windows. 
% With octave 3.8.1, the pdf does not look very nice --- it is full of "scratches". 
%print(hfig, '-depsc', 'rosenbrock');
%system('epstopdf rosenbrock.eps');