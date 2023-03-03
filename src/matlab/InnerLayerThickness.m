% Boundary layer flow
% Copyright © 2023 Ole Lindberg

function li = InnerLayerThickness(L,kappa,z0,tol,iterMax)

li = L;
iter = 0;
while(true)
    liNew = L*(2*kappa^2)/log(li/z0);
    iter = iter + 1;
    delta = abs(liNew-li);
    li = liNew;
    if (delta<tol || iter>iterMax)
        break
    end
end


