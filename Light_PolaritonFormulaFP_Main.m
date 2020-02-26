%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 25Feb2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the anti-crossing between a FabryPerot cavity mode and a 
% dipole using the Reflection of a FabryPerot formula. It sweeps over the cavity
% length in order to compute the anti-crossing

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1=3;                             %% index material 1 for Bragg mirror
n2=3.6;                           %% index material 2 for Bragg mirror
n33=3.6;                           %% index material for the FabryPerot cavity
Nperiod = 20;                     %% amount of Bragg mirror periods

N=1e23;                           %% electron density 
z0=2e-9;                          %% dipole of the transition (m)
lambda0 =950e-9;                  %% Cavity Central wavelength (m)
DeltaE=0.001;                     %% half broadening of the transition Energy eV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=(850:0.2:1050)*1e-9;         %% wavelength vector (m)
l3 = linspace(120,140,500)*1e-9;  %% cavity length vector (m)
l33=lambda0/(2*abs(n33));
l3=sort([l3 l33]);

l1=lambda0/(4*abs(n1));   % thickness at lambda/4
l2=lambda0/(4*abs(n2));   % thickness at lambda/4

w=2*pi*c./lambda;                 %% transformation of lambda in pulsation
w0=2*pi*c/lambda0;                %% transformation of lambda in pulsation
G0=2*pi*DeltaE*e/h;               %% broadening of the transition

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

Ki=(n33^2-1) + Cst*Lorentz1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Epsi=1+Ki;
R=real(Epsi);
I=imag(Epsi);

nn = (1/sqrt(2))*sqrt(R+sqrt(R.^2+I.^2));
kk = (1/sqrt(2))*sqrt(-R+sqrt(R.^2+I.^2));

n3=nn+1i*kk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LLambda,L3] = meshgrid(lambda,l3);
[N3,L3] = meshgrid(n3,l3);

alpha3=4*pi*imag(N3)./LLambda;

r12=(n1-n2)/(n1+n2);
r21=(n2-n1)/(n1+n2);

t12=2*n1/(n1+n2);
t21=2*n2/(n1+n2);

D1(1,1,:)= exp(+1i*2*pi*n1*l1./lambda);       %% take care on the sign here
D1(2,2,:)= exp(-1i*2*pi*n1*l1./lambda);       %% take care on the sign here

D2(1,1,:)= exp(+1i*2*pi*n2*l2./lambda);       %% take care on the sign here
D2(2,2,:)= exp(-1i*2*pi*n2*l2./lambda);       %% take care on the sign here

P1=(1/t12)*[1 r12 ; r12 1];
P2=(1/t21)*[1 r21 ; r21 1];

for j=1:length(lambda)
  S(:,:,j)=D2(:,:,j)*P2*D1(:,:,j)*P1;
end

for j=1:length(lambda)
  SN(:,:,j)=S(:,:,j)^Nperiod;
end

for j=1:length(lambda)
  Rb(j)=(abs(SN(1,2,j)/SN(2,2,j)))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TT=1-Rb;
theta=pi;
phi=2*pi*N3.*L3./LLambda;
delta=2*(phi-theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see the Chapther of Vincenzo Savona for the formula in the book:
% Confined Photon Systems
% Fundamentals and Applications Lectures from the Summerschool Held in Cargèse, Corsica, 3–15 August 1998
% Linear Optical Properties of Semiconductor Microcavities with Embedded Quantum Wells, Vincenzo Savona, pages 173-242
% 3) The Fabry-Perot resonator p184
% https://link.springer.com/chapter/10.1007/BFb0104383
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emmanuel Rosencher, Optoelectronic
% Complement to Chapter 9
% 9D) Fabry-Perot cavities and Bragg reflectors, page 437
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tfp = TT ./ ( 1 -  Rb .* exp(2i*(phi-theta)) ) ;
Tfp = (abs(tfp)).^2 .* exp(-alpha3.*L3);

rfp = - ( (sqrt(Rb)-sqrt(Rb).*exp(2i*(phi-theta))) ./ (1-Rb.*exp(2i*(phi-theta))));
Rfp = (abs(rfp)).^2;

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

pcolor(lambda*1e9,l3*1e9,Tfp)
plot(lambda*1e9,Tfp)
 
caxis([0 1.01])
xlim([lambda(1) lambda(end)]*1e9)
ylim([l3(1) l3(end)]*1e9)
shading flat
xlabel('Wavelength (nm)')
ylabel('Cavity length (nm)')
zlabel('Transmission')
title('Polariton anti-crossing')
colormap(jet)
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%