clear all;
close all;
clc;
format long;

a = 0;
b = 1;
N = 10;
h = (b-a)/N;
x = a:h:b;


A = zeros(N);
B = zeros(N,1);

A(1,1) = 1;
A(N+1,N+1) = 1;


for k = 2:N
    A(k,k) = trap(k,k);
    A(k,k+1) = trap(k,k+1);
    A(k,k-1) = trap(k,k-1);
end

B(1) = 0;
B(N+1) = 0;
for k = 2:N
    B(k) = trapf(k);
end

X = A\B;

exact = zeros(N+1,1);
for l = 1:N+1
    exact(l) = -x(l)*(x(l)-1)*0.5;
end
[X exact]
plot(x,X,x,exact);
