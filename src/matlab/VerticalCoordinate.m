function Z = VerticalCoordinate(z0,li,lo,eta)

a = 0;
b = 2*(li+lo); 

param.z0 = z0;
param.li = li;
param.lo = lo;

Z = 0*eta;
for i=1:length(eta)
    Z(i) = VerticalCoordinateBisectionIteration(param,eta(i),a,b,eps,100);
end

