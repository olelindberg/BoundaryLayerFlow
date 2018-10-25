function c = VerticalCoordinateBisectionIteration(param,eta,a,b,tol,iterMax)

z0 = param.z0;
li = param.li;
lo = param.lo;

a = 0;
b = 2*(li+lo); 
iter = 0;
while(true)

    c     = (a+b)/2;
    f     = log((c+z0)/z0)/log((li+z0)/z0) + c/lo - eta;

    delta = abs(a-c);
    iter  = iter + 1;
    if (delta<tol || iter>iterMax)
        break
    end

    if (f>0)
        b=c;
    else
        a=c;
    end

end

