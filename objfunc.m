clear all;
close all;

s = 1.0;

X1 = [0 : s : 10];
X2 = [-30 : s : 10];
[X1,X2] = meshgrid(X1,X2);
Z = exp((X2-X1)/2) + (5/8)*(X1-X2);

% Use log scale of Z for a batter usage of the color spectrum  
minZ = min(Z(:));  
maxZ = max(Z(:));
C = minZ + (maxZ-minZ).*log(1+Z-minZ)./log(1+maxZ-minZ); 

hfig = figure(1);
surf(X1, X2, Z, C, 'EdgeColor', 'none', 'LineStyle', 'none');
colormap = jet;

axis([0, 10, -30, 10, 0, maxZ]);
xlabel('x1', 'fontsize', 18);
ylabel('x2', 'fontsize', 18);
zlabel('f', 'fontsize', 18);