clear all
close all

showTopography              = 1;
showTopographyFourierSpace  = 0;
showVerticalCoordinate      = 0;
showReferenceVelocity       = 0;
showFourierSpaceSolutions   = 0;

isInvicid       = 0;
waveNumberMax   = 100000;

rho   = 1.25;
z0    = 1.0;
kappa = 0.4;

%------------------------------------------------------------------------------
% Horizontal coordinate:
%------------------------------------------------------------------------------
Lx = 1000;
Nx = 64;
dx = Lx/Nx;
x  = (dx/2:dx:Lx-dx/2)';

%------------------------------------------------------------------------------
% Topography:
%------------------------------------------------------------------------------
f1 = Topography(1,Lx,x);

%------------------------------------------------------------------------------
% Vertical coordinate computational space:
%------------------------------------------------------------------------------
Nz      = 256;
etaMax  = 1;
deta    = etaMax/Nz;
etaf    = (0:deta:etaMax-deta)';  
etac    = (deta/2:deta:etaMax-deta/2)';  

f1FFT = fftshift(fft(f1));
kFFT  = 2*pi/Lx*(-Nx/2:Nx/2-1)';

u1FFT = zeros(Nx,Nz);
w1FFT = zeros(Nx,Nz);

for kId = 1:min(length(kFFT),waveNumberMax)

    k1 = kFFT(kId);
    
    disp(['solving for wave number ' num2str(k1)])
    fflush(stdout); 

    %------------------------------------------------------------------------------
    % Vertical grid:
    %------------------------------------------------------------------------------
    [Zf,Zc] = VerticalGrid(Lx,kappa,z0,k1,etaf,etac);

    %------------------------------------------------------------------------------
    % Friction velocity:
    %------------------------------------------------------------------------------
    z10 = z0+10;
    u10 = 10;  % wind velocity in z=10m
    us  = u10*kappa/(log((z10+z0)/z0));

    %------------------------------------------------------------------------------
    % Reference velocity:
    %------------------------------------------------------------------------------
    [u0,u0z] = LogarithmicVelocityProfile(us,kappa,z0,Zf);
    if isInvicid
        u0  = 0*Zf+u10;
        u0z = 0*Zf;
    end

    %------------------------------------------------------------------------------
    % First order Fourier transformed Reynolds equations (eq. 6 page 276):
    %------------------------------------------------------------------------------

    rhoInv = 1/rho;

    dZc = Zc(2:end) - Zc(1:end-1);
    dZf = Zf(2:end) - Zf(1:end-1);

    Dzc = spdiags([[-1./dZc;0],[0;1./dZc]],[-1 0],Nz,Nz);
    Dzf = spdiags([[-1./dZf;0],[0;1./dZf]],[ 0 1],Nz,Nz);

    ZERO = 0*speye(Nz,Nz);

    % Momentum conservation in x-direction:
    Eq11 = spdiags(1i*k1*u0    ,0,Nz,Nz);
    Eq12 = spdiags(u0z         ,0,Nz,Nz);
    Eq13 = spdiags(rhoInv*1i*k1*ones(Nz,1),0,Nz,Nz);
    if isInvicid
        Eq14 = ZERO;
    else
        Eq14 = -Dzc;
    end

    % Boundary conditions on u1:
    Eq11(1,:) = 0;
    Eq12(1,:) = 0;
    Eq13(1,:) = 0;
    Eq14(1,:) = 0;
    Eq11(1,1) = 1;
    % Assemble:
    Eq1  = [Eq11,Eq12,Eq13,Eq14];
    rhs1 = zeros(Nz,1);

    % Momentum conservation in z-direction:
    Eq21 = ZERO;
    Eq22 = spdiags(1i*k1*u0    ,0,Nz,Nz);
    Eq23 = 1/rho*Dzc;
    Eq24 = ZERO;
    % Boundary conditions on w1:
    Eq21(1,:) = 0;
    Eq22(1,:) = 0;
    Eq23(1,:) = 0;
    Eq24(1,:) = 0;
    Eq22(1,1) = 1;
    % Assemble:
    Eq2  = [Eq21,Eq22,Eq23,Eq24];
    rhs2 = (k1*u0).^2*f1FFT(kId );
    rhs2(1) = 0;

    % Mass conservation:
    Eq31 = spdiags(1i*k1*ones(Nz,1),0,Nz,Nz);
    Eq32 = Dzf;
    Eq33 = ZERO;
    Eq34 = ZERO;
    % Boundary conditions on p1:
    Eq31(end,:) = 0;
    Eq32(end,:) = 0;
    Eq33(end,:) = 0;
    Eq34(end,:) = 0;
    
    alphak = sqrt(k1*k1);
    uop    = k1*u0(end);
    Eq32(end,end) = -1i*uop/alphak;
    Eq33(end,end) = 1;
    % Assemble:
    Eq3  = [Eq31,Eq32,Eq33,Eq34];
    rhs3 = zeros(Nz,1);
    rhs3(end) = - uop^2/alphak*f1FFT(kId );

    % Mixing-length closure:
    Eq41 = spdiags(2*kappa*(Zc-z0)*us,0,Nz,Nz)*Dzf;
    Eq42 = ZERO;
    Eq43 = ZERO;
    Eq44 = -speye(Nz,Nz);
    % Boundary conditions on tau1:
    Eq41(end,:) = 0;
    Eq42(end,:) = 0;
    Eq43(end,:) = 0;
    Eq44(end,:) = 0;
    Eq44(end,end) = 1;
    % Assemble:
    Eq4  = [Eq41,Eq42,Eq43,Eq44];
    rhs4 = zeros(Nz,1);

    % Assemble matrix equations:
    A   = [Eq1 ;Eq2 ;Eq3 ;Eq4 ];
    rhs = [rhs1;rhs2;rhs3;rhs4];

    % Solve:
    sol = A\rhs;
    sol = reshape(sol,Nz,4);

    % Extract solution:
    u1   = sol(:,1); 
    w1   = sol(:,2);
    p1   = sol(:,3);
    tau1 = sol(:,4);

    u1FFT(kId,:) = u1;
    w1FFT(kId,:) = w1;

    if showFourierSpaceSolutions
        figure
        subplot(2,2,1)
        plot(u1,Zf)
        grid on
        subplot(2,2,2)
        plot(w1,Zf)
        grid on
        subplot(2,2,3)
        plot(p1,Zc)
        grid on
        subplot(2,2,4)
        plot(tau1,Zc)
        grid on
    end
end

u1All = 0*u1FFT;
w1All = 0*w1FFT;
for i = 1:Nz
    u1All(:,i) = ifft(fftshift(u1FFT(:,i)));
    w1All(:,i) = ifft(fftshift(w1FFT(:,i)));
end
u1All = real(u1All);
w1All = real(w1All);

%------------------------------------------------------------------------------
% Plots:
%------------------------------------------------------------------------------
figure
for i = 1:Nx
    plot(x(i)+u1All(i,:),f1(i)+Zf-z0)
    hold on
end
if showTopography
    plot(x,f1,'k-','LineWidth',2)
    grid on
    xlabel('x')
    ylabel('z')
end

if showTopographyFourierSpace
    figure
    plot(kFFT,real(f1FFT),'r-','LineWidth',2)
    hold on
    plot(kFFT,imag(f1FFT),'g-','LineWidth',2)
    grid on
    xlabel('x')
    ylabel('z')
end


if showVerticalCoordinate
    figure
    plot(etaf,Zf,'r-','LineWidth',2)
    hold on
    plot(etac,Zc,'g-','LineWidth',2)
    plot(etaf,Zf,'b.')
    plot(etac,Zc,'c.')
    title('Coordinate transform')
    xlabel('\eta')
    ylabel('Z')
    grid on
end

if showReferenceVelocity
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
end


