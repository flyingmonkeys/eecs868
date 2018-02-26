clear all;
close all;

a_max = 10.0;

c1 = 0.6;
c2 = 0.7;
x = [1 3]';
d = [1/sqrt(2) -1/sqrt(2)]';
a = [0:0.1:a_max];

for k=1:length(a)
	phi(k) = f(x + a(k)*d); 		% compute phi
	arm(k) = f(x) + c1*a(k)*phiprime(0,d,x); 	% compute Armijo line
	curv(k) = phiprime(a(k),d,x);		% compute curvature line
end


figure(1)
plot(a,phi);
hold on;
plot(a,arm);

figure(2)
plot(a,c2*abs(phiprime(0,d,x))*ones(length(a)));
hold on;
plot(a,abs(curv));