clear all;
clc;
%Declaring the 'a' Value of the ray
a = 2;
%Finding the Critical Radius
R = (18*18/12400)^(-1/6);
%Declaring the function theta
%theta = @(r) a./(r*((1-38.27*r.^(-4))-a^2));
theta = @(r) a./(((r.^4)-38.27*(r.^-2)-a^2*(r.^2)));
fun = @(x) x.^4 - 38.27*x.^-2 - a^2*x.^2;
roots = fzero(fun,1);
%Calculating the critical angle theta_critical
theta_critical = integral(theta,roots,Inf);
%Defining theta matrix
theta_val = zeros(1,10000);
%Intial Value of R
Ri = sqrt(a^2 + 5^2);
R_val = linspace(Ri,1.94d,10000);
for i = 1:10000
    theta_val(i) = integral(theta,R_val(i),inf);
end
polar(theta_val,R_val);

