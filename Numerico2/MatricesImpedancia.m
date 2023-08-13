clear all; close all; clc; 

a           =   5/(2*pi); % Radio del círculo en longitud de onda
phi_i       =   0; % Angulo de incidencia en grados
r_=5; N=100; eta0=120*pi; % Rango del grafico en longitudes de onda, numero de intervalos y constante

[Z_TM,I_TM,RCS_TM,Z_TE,I_TE,RCS_TE,phi,E,Es,Ei,H,Hs,Hi]=RCSedit(a,N,phi_i,r_);


figure()
imagesc(abs(Z_TE));
colorbar; 
title('Matriz de impedancia para TE','Interpret','Latex','FontSize',14);
xlabel('Columnas'); ylabel('Filas');

figure()
imagesc(abs(Z_TM));
colorbar; 
title('Matriz de impedancia para TM','Interpret','Latex','FontSize',14);
xlabel('Columnas'); ylabel('Filas');