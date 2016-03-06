function V = trapf(i)

sum = 0;
a = 0;
b = 1;
n = 1000;
dx = (b-a)/n;
pos = a:dx:b;


for k  = 1:n
    %sum = sum + (pos(k+1)-pos(k))*0.5*((f(pos(k))*phi(pos(k),i)) + (f(pos(k+1))*phi(pos(k+1),i)));
    sum = sum + dx*((phip(pos(k),i)*f(pos(k))));
end

V =sum;