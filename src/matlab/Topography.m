% Boundary layer flow
% Copyright © 2023 Ole Lindberg

function f1 = Topography(topographyId,Lx,x)

TOPOGRAPHY.COSINE   = 0;
TOPOGRAPHY.GAUSSIAN = 1;
TOPOGRAPHY.PLATEAU  = 2;

if topographyId == TOPOGRAPHY.COSINE
    A1 = 50;
    f1 = A1*cos(2*pi/Lx*x);
elseif topographyId == TOPOGRAPHY.GAUSSIAN
    x0    = 500;
    sigma = 100;
    H     = 100;
    f1    = GaussianFunction(H,x0,sigma,x);
elseif topographyId == TOPOGRAPHY.PLATEAU
    Nx      = length(x);
    sigma   = 100;
    H       = 100;
    f1      = 0*x;
    id     = 1:Nx/4;
    f1(id) = GaussianFunction(H,Lx/4,sigma,x(id));
    id     = Nx/4+1:Nx*3/4;
    f1(id) = H;
    id     = Nx*3/4+1:Nx;
    f1(id) = GaussianFunction(H,Lx*3/4,sigma,x(id));
else
    error('No such topography')
end

