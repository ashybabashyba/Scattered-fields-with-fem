clear all; close all; clc; 

a           =   5/(2*pi); % Radio del círculo en longitud de onda
phi_i       =   0; % Angulo de incidencia en grados
r_=5; N=100; eta0=120*pi; % Rango del grafico en longitudes de onda, numero de intervalos y constante

[Z_TM,I_TM,RCS_TM,Z_TE,I_TE,RCS_TE,phi,E,Es,Ei,H,Hs,Hi]=RCSedit(a,N,phi_i,r_);
[x,y]   =   meshgrid(linspace(-r_,r_,N),linspace(-r_,r_,N));

figure()
pcolor(x,y,abs(E))
grid on;
xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$y/\lambda$','Interpret','Latex','FontSize',14)
title('$|E|$ para TM','Interpret','Latex','FontSize',14)
colorbar ; axis equal; colormap jet;
set(gca,'TickLabel','Latex','FontSize',15)
set(colorbar,'TickLabelInterpret','Latex','FontSize',14)
axis([-1 +1 -1 +1]*r_); caxis([0 2]);


figure()
pcolor(x,y,abs(H)*eta0)
xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$y/\lambda$','Interpret','Latex','FontSize',14)
title('$\eta|H|$ para TE','Interpret','Latex','FontSize',14)
colorbar; axis equal; colormap jet;
set(gca,'TickLabel','Latex','FontSize',15)
set(colorbar,'TickLabelInterpret','Latex','FontSize',14)
axis([-1 +1 -1 +1]*r_); caxis([0 2])

