function V = phi(c,i)

a = 0;
b = 1;
N =10;
h = (b-a)/N;
x = a:h:b;

for j=1:N+1
    if (i==1)
        V = (c - a)/h;
    elseif (i==N+1)
        V = (b - c)/h;
    elseif ( c>= x(i-1) && c<=x(i))
        V = (c - x(i-1))/(x(i) - x(i-1));
    elseif ( c>= x(i) && c<=x(i+1) )
        V = (x(i+1) - c)/(x(i+1) - x(i));
    else
        V = 0;
    end
end
