clear all; close all; clc;
a           =   5/(2*pi); % Radio del círculo en longitud de onda
phi_i       =   180; % Angulo de incidencia en grados
colors=jet(90); r_=5; % Rango del grafico en longitudes de onda, numero de intervalos y constante

for N=10:20:90
    figure(1);
    hold on
    [Z_TM,I_TM,RCS_TM,Z_TE,I_TE,RCS_TE,phi,E,Es,Ei,H,Hs,Hi]=RCSedit(a,N,phi_i,r_);
    stairs(phi*360/pi,abs(I_TM),'LineWidth',1, 'color',colors(N,:),'DisplayName',strcat(num2str(N), ' celdas o segmentos.'))
    axis([0 180 0 2]); grid on; 
    hold off;
end
legend('show'); xlabel('Ángulo \phi'); ylabel('Magnitud de la densidad de corriente |J|')
