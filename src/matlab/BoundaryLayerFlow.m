clear all
close all

rho   = 1.25;
z0    = 0.03;
kappa = 0.4;

%------------------------------------------------------------------------------
% Horizontal coordinate:
%------------------------------------------------------------------------------
Lx = 1000;
Nx = 32;
dx = Lx/Nx;
x  = 0:dx:Lx;

%------------------------------------------------------------------------------
% Topography:
%------------------------------------------------------------------------------
k1 = 2*pi/Lx;
A1 = 50;
f1 = A1*cos(k1*x);

%------------------------------------------------------------------------------
% Inner layer thickness:
%------------------------------------------------------------------------------
li = InnerLayerThickness(Lx,kappa,z0,1e-10,1000);

%------------------------------------------------------------------------------
% Outer layer thickness:
%------------------------------------------------------------------------------
lo = OuterLayerThickness(k1);

%------------------------------------------------------------------------------
% Vertical coordinate:
%------------------------------------------------------------------------------
Nz      = 32;
etaMax  = 3;
deta    = etaMax/Nz;
etaf    = 0:deta:etaMax;  
etac    = deta/2:deta:etaMax-deta/2;  

Zf = VerticalCoordinate(z0,li,lo,etaf);
Zc = VerticalCoordinate(z0,li,lo,etac);

%------------------------------------------------------------------------------
% Friction velocity:
%------------------------------------------------------------------------------
z10 = 10;
u10 = 10;  % wind velocity in z=10m
us = u10*kappa/(log((z10+z0)/z0));

%------------------------------------------------------------------------------
% Reference velocity:
%------------------------------------------------------------------------------
[u0,u0z] = LogarithmicVelocityProfile(us,kappa,z0,Zf);

%------------------------------------------------------------------------------
% First order Fourier transformed Reynolds equations (eq. 6 page 276):
%------------------------------------------------------------------------------

% Boundary condition:
k = 1;

k   = 2;
dZc = Zc(k)   - Zc(k-1);
dZf = Zf(k+1) - Zf(k);

%      u1         , w1         , p1           , tau1    

A   = [0          , 0          , 0.5*1i*k1/rho, 1/dZc  ;
       0          , 0          , 1/dZc        , 0      ;
       0          , 0          , 0            , 0      ];

B   = [1i*k1*u0(k), u0z(k)     , 0.5*1i*k1/rho, -1/dZc ;
       0          , 1i*k1*u0(k), -1/dZc       , 0      ;
       0.5*1i*k1  , -1/dZf     , 0            , 0      ];

C   = [0          , 0          , 0            , 0      ;
       0          , 0          , 0            , 0      ;
       0.5*1i*k1  , 1/dZf      , 0            , 0      ];





%------------------------------------------------------------------------------
% Plots:
%------------------------------------------------------------------------------
figure
plot(x,f1,'k-','LineWidth',2)
grid on
axis equal
xlabel('x')
ylabel('z')

figure
plot(etaf,Zf,'r-','LineWidth',2)
hold on
%plot(etac,Zc,'g--','LineWidth',2)
plot(etaf,Zf,'b.')
plot(etac,Zc,'c.')
title('Coordinate transform')
xlabel('\eta')
ylabel('Z')
grid on

figure
plot([0 12],[li li],'k-','LineWidth',2)
hold on
plot([0 12],[li+lo li+lo],'k-','LineWidth',2)
plot(u0,Zf,'r-','LineWidth',2)
plot(u0,Zf,'b.')
title('Logarithmic velocity profile')
xlabel('u_0')
ylabel('Z')
grid on



