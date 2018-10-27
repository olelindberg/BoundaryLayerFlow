function [Zf,Zc] = VerticalGrid(Lx,kappa,z0,k1,etaf,etac)

verticalGridId = 0;

if verticalGridId==0
    zMax = 0.10*Lx;
    Zf = etaf*zMax;
    Zc = etac*zMax;
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
