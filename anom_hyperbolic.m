function [F] = anom_hyperbolic(M,e)

format long
syms x;
f=e*sinh(x)-x-M; %Enter the Function here
g=diff(f); %The Derivative of the Function

n=11;

epsilon = 10^-(n);

x0 = M;
for i=1:100
     f0=vpa(subs(f,x,x0)); %Calculating the value of function at x0
     f0_der=vpa(subs(g,x,x0)); %Calculating the value of function derivative at x0
    y=x0-f0/f0_der; % The Formula
    err=abs(y-x0);
    
    if err<epsilon %checking the amount of error at each iteration
        break
    end
    x0=y;
end
y = y - rem(y,10^-n); %Displaying upto required decimal places

F=double(y);
end