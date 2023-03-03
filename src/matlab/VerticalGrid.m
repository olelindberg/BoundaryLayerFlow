% Boundary layer flow
% Copyright © 2023 Ole Lindberg

function [Zf,Zc] = VerticalGrid(Lx,kappa,z0,k1,etaf,etac)

verticalGridId = 0;

if verticalGridId==0
    zMax = 0.20*Lx;
    Zf = z0+etaf*(zMax-z0);
    Zc = z0+etac*(zMax-z0);
elseif verticalGridId==1
    %------------------------------------------------------------------------------
    % Inner layer thickness:
    %------------------------------------------------------------------------------
    li = InnerLayerThickness(Lx,kappa,z0,1e-10,1000);

    %------------------------------------------------------------------------------
    % Outer layer thickness:
    %------------------------------------------------------------------------------
    lo = OuterLayerThickness(k1);

    %------------------------------------------------------------------------------
    % Vertical coordinate physical space:
    %------------------------------------------------------------------------------
    Zf = VerticalCoordinate(z0,li,lo,etaf);
    Zc = VerticalCoordinate(z0,li,lo,etac);
end
