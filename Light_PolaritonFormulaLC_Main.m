%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 26Feb2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the anti-crossing between a LC cavity mode and a dipole 
% using the impedance line formula. It sweeps over the inductance length in 
% order to compute the anti-crossing

% "Interaction between meta-materials and shallow donors in bulk GaN at THz frequency"
% Laurent Nevou, Etienne Giraud, Fabrizio Castellano, Nicolas Grandjean, and Jerome Faist 
% Optics Express, 22, 3, p3199 (2014)

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 6.62606896e-34;           %% Planck constant (J.s)
hbar  = h/(2*pi);
e     = 1.602176487e-19;          %% Electron charge (Coulomb)
c     = 2.99792458e8;             %% Speed of light (m/s)
Epsi0 = 8.854187817620e-12;       %% Vaccum dielectric constant (F/m)
mu0   = 1/(Epsi0*c^2);            %% Vaccum permeabiliy (A/m)
Z0    = mu0*c;                    %% Vaccum impedance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LC parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wd = 0.9e-6;       %% width of the inductance (meter)
t  = 70e-9;        %% thickness of the inductor (meter)
d  = 2e-6;         %% distance between the 2 plates of the capacitor
l1 = 2e-6;         %% dimension of the capcitor plate
l2 =0.84e-6;       %% dimension of the capcitor plate
A  = l1*l2;        %% Area of the capacitor plate
nc = 2.5;          %% Substrate optical refraction index at high frequency
Er = nc^2;
ZG = Z0/nc;        %% substrate impedance
duty=0.53;         %% duty cycle of resonator density over the substrate
ro=2.44E-8;        %% metal resistivity (Ohm.meter) (ro=2.44E-8 Ohm.m for Au and 1.68E-8 Ohm.m for Copper)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dipole parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=1e23;                           %% electron density 
z0=2e-9;                          %% dipole of the transition (m)
lambda0 =100e-6;                  %% Cavity Central wavelength (m)
DeltaE=0.0005;                    %% half broadening of the transition Energy eV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=(50:0.1:150)*1e-6;         %% wavelength vector (m)

w = 2*pi*c./lambda;               %% transformation of lambda in pulsation
w0= 2*pi*c/lambda0;               %% transformation of lambda in pulsation
G0= 2*pi*DeltaE*e/h;              %% broadening of the transition
f = w/2/pi;                       %% transformation of the pulsation in frequency

l = linspace(0,120e-6,500);       %% Variable definition length of the inductance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Lorentzian dipole build %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A. Gabbay, J. Reno, J. R. Wendt, A. Gin, M. C. Wanke, M. B. Sinclair, E. Shaner, and I. Brener,
% “Interaction between metamaterial resonators and intersubband transitions in semiconductor quantum wells,” 
% Appl. Phys. Lett. 98(20), 203103 (2011). 

% J. Faist, "Optical properties of semiconductor"
% chap3: "Light-matter interaction"
% 3.5.2 A polarization field

Cst=N*e^2*z0^2/(Epsi0*hbar);

Lorentz1 = w0 ./ ( w0^2-w.^2 - 1i*G0*w );
%Lorentz2 = 1/2 * ( w0-w+1i*G0/2 ) ./ ( (w0-w).^2 + (G0/2)^2 );
%Lorentz3 = 1/2 * (1./(w+w0+1i*G0/2) - 1./(w-w0+1i*G0/2) );

Ki=(nc^2-1) + Cst*Lorentz1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Epsi=1+Ki;
R=real(Epsi);
I=imag(Epsi);

nn = (1/sqrt(2))*sqrt(R+sqrt(R.^2+I.^2));
kk = (1/sqrt(2))*sqrt(-R+sqrt(R.^2+I.^2));

nc=nn+1i*kk;
%alpha=2*kk.*w/c;       % Absorbance (m-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Capacitance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://en.wikipedia.org/wiki/Capacitance

c = (Epsi0 * Epsi * A) / d  ;     %% Capacitor formula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FF,L] = meshgrid(f,l);
[CC,L] = meshgrid(c,l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inductance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Signal Integrity and Radiated Emission of High?Speed Digital Systems"
% Spartaco Caniggia Francescaromana Maradei
% DOI:10.1002/9780470772874
% Appendix A: Formulae for Partial Inductance Calculation (Pages: 481-486)
% https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470772874.app1

LL = mu0/2/pi * L .*( log(2*L/(wd+t)) + 1/2 + 2/9*(wd+t)./L ) ; %% inductance (Henry)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Resistance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = ro * L/(wd*t);                                 %% stripe resistivity (Ohm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Impedance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z_RLC = RR - 1i*2*pi*FF.*LL - 1./(1i*2*pi*FF.*CC);  %% RLC complex impedance

B=Z0*( 1/ZG + 1./Z_RLC );                           %% Normalization constante
T0=4*(Z0*ZG)/((Z0+ZG)^2);                           %% Substrate transmission

trans1=duty*4*(Z0/ZG)*(1./(1+B)).*conj(1./(1+B))  +  (1-duty)*T0;
trans1=trans1/T0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X0fig=-3500; Y0fig=100;
X0fig=100; Y0fig=100;
Wfig=1500;Hfig=900;

figure('Name','Results','position',[X0fig Y0fig Wfig Hfig])
FS=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on; box on;

pcolor(lambda*1e6,l*1e6,abs(trans1))

shading flat
colormap(jet)
colorbar

%caxis([0 1])
xlim([lambda(1) lambda(end)]*1e6)
ylim([l(1) l(end)]*1e6)
xlabel('Wavelength (um)')
ylabel('Inductance length (µm)')
zlabel('Transmission')
title('Polariton anti-crossing')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%