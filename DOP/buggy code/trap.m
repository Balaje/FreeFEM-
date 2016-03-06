function V = trap(i,j)

sum = 0;
a = 0;
b = 1;
n = 1000;
dx = (b-a)/n;
pos = a:dx:b;

for k = 1:n
    %sum = sum + (pos(k+1)-pos(k))*0.5*((phip(pos(k),i)*phip(pos(k),j)) + (phip(pos(k+1),i)*phip(pos(k+1),j)));
    sum = sum + dx*((phip(pos(k),i)*phip(pos(k),j)));
end

V = sum;