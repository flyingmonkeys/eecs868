% Wolfe.m
% Plotting script to help visualize the Wolfe conditions (Armijo and curvature) vs. f

clear all;
close all;

a_max = 4;

c1 = 0.4;
c2 = 0.7;
x = [0.44059 .19011]';
d = [.456 .890]';
a = [0:0.05:a_max];

for k=1:length(a)
	phi(k) = f(x + a(k)*d); 		% compute phi
	arm(k) = f(x) + c1*a(k)*phiprime(0,d,x); 	% compute Armijo line
	curv(k) = phiprime(a(k),d,x);		% compute curvature line
end


figure(1)
plot(a,phi);
hold on;
plot(a,arm,'--');
title('\phi vs. \alpha with superimposed Armijo constraint (c_1 = 0.01)');
xlabel('\alpha');
ylabel('\phi');
grid on;
legend({'\phi','Armijo constraint'});
%axis([0 4 -1.5 3]);

figure(2)
plot(a,abs(curv));
hold on;
plot(a,c2*abs(phiprime(0,d,x))*ones(length(a)),'--');
title('|\phi''| vs. \alpha with superimposed curvature constraint (c_2 = 0.10)');
xlabel('\alpha');
ylabel('|\phi''|');
grid on;
legend({'|\phi''|','Curv. constraint'});
%axis([0 4 -1.2 1.2]);