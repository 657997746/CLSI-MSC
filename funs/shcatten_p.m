function x = shcatten_p(a, r, p)

x = 0;
if p == 1  
    if a > r
        x = a-r;
    end;
elseif p == 0
    if a > sqrt(2*r)
        x = a;
    end;
else
    v = (r*p*(1-p))^(1/(2-p))+eps;
v1 = v+r*p*v^(p-1);
ob0 = 0.5*a^2;
if a > v1
    x = a;
    for i = 1:10
        f = (x-a) + r*p*x^(p-1);
        g = 1-r*p*(1-p)*x^(p-2);
        x = x-f/g;
    end;
    ob1 = 0.5*(x-a)^2 + r*x^p;
    x_can = [0,x];
    [temp,idx] = min([ob0,ob1]);
    x = x_can(idx);
end;
end;