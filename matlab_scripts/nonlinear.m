clear all;
clc;
%To find and plot R(theta)
%The two boundary conditions are : 
%1) R(tan^-1(c/5R)) = sqrt(c^2 + 25R^2)
%2) R'(tan^-1(c/5R)) = 0

c = 0.5;                  %'c' Parameter
lb = atan(c/10);         %lb on theta
ub = pi;              %ub on theta (Till the Singularity)

delta = 0.001;          %grid width (Accuracy of the problem

theta = lb:delta:ub;    %theta matrix

l = length(theta);      %size of the grid

r = zeros(1,l);           %initializing the R vector

%Incorporating the Boundary Conditions
r(1) = sqrt(c^2 + 100);
r(2) = r(1)*(1-delta*cos(lb));
%r(2) = r(1);

%Defining the Difference Equation for the given Problem
for k = 2:l-1
    %Defining the non-linear function f(r)
    f = (((7400)*(1+r(k)^2))/(324*(r(k)^6)-12400));
    %Defining the quadratic coefficients of x(k+1)
    a = r(k) - delta;
    b = (-3*(r(k)^2) + r(k)*r(k-1) + 2*delta*r(k));
    c = (2*(r(k)^3) - (r(k)^2)*r(k-1) - (delta * (r(k)^2)) - (delta^3)*f);
    r(k+1) = (((-1*b)-sqrt(b^2-4*a*c))/(2*a));
end

%Plotting the Polar Plot of R(theta)
polar(theta,r);
