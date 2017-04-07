clear all;
clc;
%Declaring the 'a' Value of the ray
a = 0.1;
%Finding the Critical Radius
R = (18*18/12400)^(-1/6);
%Declaring the function theta
%theta = @(r) a./(r*((1-38.27*r.^(-4))-a^2));
syms r
theta =  a./(((r.^4)-38.27*(r.^-2)-a^2*(r.^2)));
%syms x;
fun =@(x) x.^4 - (38.27)*x.^-2 - a^2*x.^2;
ezplot(theta);
roots = fzero(fun,2);
%plot(fun);
%Calculating the critical angle theta_critical
%theta_critical = integral(theta,roots,Inf);
%Defining theta matrix
%theta_val = zeros(1,10000);
%Intial Value of R
%Ri = sqrt(a^2 + 5^2);
%R_val = linspace(Ri,roots,10000);
%for i = 1:10000
%    theta_val(i) = integral(theta,R_val(i),inf);
%end
%polar(theta_val,R_val);

