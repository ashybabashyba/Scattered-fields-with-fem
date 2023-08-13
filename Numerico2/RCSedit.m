function[Z_TM,I_TM,RCS_TM,Z_TE,I_TE,RCS_TE,phi,E,Es,Ei,H,Hs,Hi]=RCSedit(a,N,phi_i,r_)
Ns          =   N;
j           =   sqrt(-1);
eta0        =   120*pi;
k           =   2*pi;
%% TM-EFIE 
phi_i       =   deg2rad(phi_i);
Data        =   CircleData(a,N);
Z           =   zeros(N,N);
V           =   zeros(N,1);
for i=1:N
    for ii=1:N
        ln        	=   Data(ii,5);
        if i==ii
            func        =   @(l) besselh(0,2,k*abs(l-ln/2))+j*(2/pi)*log(l/2);
            Z(i,ii)     =   (k*eta0/4)*(Quad(@(l)func(l),0,ln)-j*(2/pi)*ln*(log(ln/2)-1));
        else
            Z(i,ii)     =   (k*eta0/4)*Quad(@(l)Integrand_TM(Data,l,i,ii,k),0,ln);    
        end
    end
end
for i=1:N
    xm          =   Data(i,1);
    ym          =   Data(i,2);
    lm          =   Data(i,5);
    phi_m     	=   Data(i,6);
    xm_         =   xm+0.5*lm*cos(phi_m);
    ym_         =   ym+0.5*lm*sin(phi_m);
    V(i,1)      =   exp(j*k*(xm_*cos(phi_i)+ym_*sin(phi_i)));
end
I           =   Z\V;
%% Determinación Campos eléctricos
Es = 0; 
[x,y]   =   meshgrid(linspace(-r_,r_,Ns),linspace(-r_,r_,Ns));
for i=1:N
    xm          =   Data(i,1);
    ym          =   Data(i,2);
    R       =   sqrt((x-xm).^2+(y-ym).^2);
    val     =   0.1;
    R       =   sqrt((x-xm).^2+(y-ym).^2).*(R>val)+sqrt((x-xm).^2+(y-ym).^2+1E-4./R).*(R<=val);
    lm          =   Data(i,5);
    Es      =   Es-(k*eta0/4)*(lm)*I(i,1)*besselh(0,2,k*R);
end
Ei      =   exp(j*k*(x*cos(phi_i)+y*sin(phi_i)));
E       =   Ei+Es;


%%
phi         =   linspace(0,pi,Ns);
sigma       =   0;
for i=1:N
    xn          =   Data(i,1);
    yn          =   Data(i,2);
    ln          =   Data(i,5);
    phi_n     	=   Data(i,6);
    xn_         =   xn+0.5*ln*cos(phi_n);
    yn_         =   yn+0.5*ln*sin(phi_n);
    sigma       =   sigma+ ln*I(i,1)*exp(j*k*(xn_*cos(phi)+yn_*sin(phi)));
end
sigma       =   (k*eta0^2/4)*abs(sigma).^2;
%%
Z_TM        =   Z;
I_TM        =   I;
RCS_TM      =   10*log10(sigma);
%% TE-MFIE
for i=1:N
    for ii=1:N
        ln         	=   Data(ii,5);
        if i==ii
            Z(i,ii)     =   -0.5;   
        else
            Z(i,ii)     =   -(j*k/4)*Quad(@(l)Integrand_TE(Data,l,i,ii,k),0,ln);   
        end
    end
end
for i=1:N
    xm          =   Data(i,1);
    ym          =   Data(i,2);
    lm          =   Data(i,5);
    phi_m     	=   Data(i,6);
    xm_         =   xm+0.5*lm*cos(phi_m);
    ym_         =   ym+0.5*lm*sin(phi_m);
    V(i,1)      =   exp(j*k*(xm_*cos(phi_i)+ym_*sin(phi_i)));
end
I           =   Z\V;


%%
phi         =   linspace(0,pi,Ns);
sigma       =   0;
for i=1:N
    xn          =   Data(i,1);
    yn          =   Data(i,2);
    ln          =   Data(i,5);
    phi_n     	=   Data(i,6);
    nx         	=   Data(i,7);
    ny         	=   Data(i,8);
    xn_         =   xn+0.5*ln*cos(phi_n);
    yn_         =   yn+0.5*ln*sin(phi_n);
    Dotn        =   (nx*cos(phi)+ny*sin(phi));
    sigma       =   sigma+ ln*I(i,1)*Dotn.*exp(j*k*(xn_*cos(phi)+yn_*sin(phi)));
end
sigma       =   (k/4)*abs(sigma).^2;
%%
Z_TE        =   Z;
I_TE        =   I;
RCS_TE      =   10*log10(sigma);


%% Determinación Campos magnéticos
Hs = 0;
[x,y]   =   meshgrid(linspace(-r_,r_,Ns),linspace(-r_,r_,Ns));
for i=1:N
    xm          =   Data(i,1);
    ym          =   Data(i,2);
    R       =   sqrt((x-xm).^2+(y-ym).^2);
    val     =   0.1;
    R       =   sqrt((x-xm).^2+(y-ym).^2).*(R>val)+sqrt((x-xm).^2+(y-ym).^2+1E-4./R).*(R<=val);
    Rx      =   (x-xm)./R;
    Ry      =   (y-ym)./R;
    lm          =   Data(i,5);
    phi_m     	=   Data(i,6);
    xm_         =   xm+0.5*lm*cos(phi_m);
    ym_         =   ym+0.5*lm*sin(phi_m);
    Cross   =   Rx*(Data(i,4)-ym_)-Ry*(Data(i,3)-xm_);
    %Hs      =   Hs-(j*k/4)*I(i,1)*(Cross/2).*besselh(1,2,k*R);
    Hs      =   Hs-(j*k*lm/4)*I_TE(i,1)*(Cross/2).*besselh(1,2,k*R);
end
Hi      =   exp(j*k*(x*cos(phi_i)+y*sin(phi_i)))/eta0;
H       =   Hi+Hs;


end
%%

function[Data]=CircleData(a,N)
%% Columnas de Data: xn yn xn+1 yn+1 ln phi_n nx ny
Data        =   zeros(N,2);
for i=1:N
    xn          =   a*cos((2*pi/N)*(i-1));
    yn          =   a*sin((2*pi/N)*(i-1));
    Data(i,1)   =   xn;
    Data(i,2)   =   yn;
end
Data        =   [ Data [Data(2:N,:);Data(1,:)] ];
ln          =   sqrt((Data(:,3)-Data(:,1)).^2+(Data(:,4)-Data(:,2)).^2);
phi_n       =   atan2((Data(:,4)-Data(:,2)),(Data(:,3)-Data(:,1)));
nx          =   (Data(:,4)-Data(:,2))./ln;
ny          =   -(Data(:,3)-Data(:,1))./ln;
Data        =   [ Data ln phi_n nx ny];
end
%%
function[I]=Integrand_TM(Data,l,m,n,k)
I           =   besselh(0,2,k*Rmn(l,Data,m,n));
end
%%
function[I]=Integrand_TE(Data,l,m,n,k)
%%
xn          =   Data(n,1);
yn          =   Data(n,2);
phi_n     	=   Data(n,6);
nx         	=   Data(n,7);
ny         	=   Data(n,8);
xm          =   Data(m,1);
ym          =   Data(m,2);
lm          =   Data(m,5);
phi_m     	=   Data(m,6);
x           =   xn+l*cos(phi_n);
y           =   yn+l*sin(phi_n);
xm_         =   xm+0.5*lm*cos(phi_m);
ym_         =   ym+0.5*lm*sin(phi_m);
Rmnx        =   (xm_-x);
Rmny        =   (ym_-y);
%%
Dotmn       =   (nx*Rmnx+ny*Rmny)./(Rmn(l,Data,m,n));       
I           =   Dotmn.*besselh(1,2,k*Rmn(l,Data,m,n));
end
%%
function[R]=Rmn(l,Data,m,n)
%%
xn          =   Data(n,1);
yn          =   Data(n,2);
phi_n     	=   Data(n,6);
xm          =   Data(m,1);
ym          =   Data(m,2);
lm          =   Data(m,5);
phi_m     	=   Data(m,6);
%%
x           =   xn+l*cos(phi_n);
y           =   yn+l*sin(phi_n);
xm_         =   xm+0.5*lm*cos(phi_m);
ym_         =   ym+0.5*lm*sin(phi_m);
%%
R           =   sqrt((xm_-x).^2+(ym_-y).^2);
end
%%
function[sigma]=PEC_RCS_TM(a,phi)
N           =   100;
k           =   2*pi;
%%
Sum         =   (besselj(0,k*a)/besselh(0,2,k*a));
for n=1:N
Sum         =   Sum+Term(n,phi,a);
end
sigma       =   (2/pi)*abs(Sum).^2;
    function[Sn]=Term(n,phi,a)
    k           =   2*pi;
    Sn          =   2*(besselj(n,k*a)/besselh(n,2,k*a))*cos(n*phi);
    end
end
%%
function[sigma]=PEC_RCS_TE(a,phi)
N           =   100;
k           =   2*pi;
%%
Sum         =   (besselj(1,k*a)/besselh(1,2,k*a));
for n=1:N
Sum         =   Sum+Term(n,phi,a);
end
sigma       =   (2/pi)*abs(Sum).^2;
function[Sn]=Term(n,phi,a)
    k           =   2*pi;
    Sn          =   2*((besselj(n-1,k*a)-(n/(k*a))*besselj(n,k*a))...
        /(besselh(n-1,2,k*a)-(n/(k*a))*besselh(n,2,k*a)))*cos(n*phi);
end
end
%%
function[I]=Quad(func,a,b)
%% 32 Levels (Gauss-Legendre rule) 
x       =   [-9.97263861849481e-01 -9.85611511545268e-01 -9.64762255587506e-01 ...
             -9.34906075937740e-01 -8.96321155766052e-01 -8.49367613732570e-01 ...
             -7.94483795967942e-01 -7.32182118740290e-01 -6.63044266930215e-01 ...
             -5.87715757240762e-01 -5.06899908932229e-01 -4.21351276130635e-01 ...
             -3.31868602282128e-01 -2.39287362252137e-01 -1.44471961582796e-01 ...
             -4.83076656877383e-02 +4.83076656877382e-02 +1.44471961582797e-01 ...
             +2.39287362252137e-01 +3.31868602282128e-01 +4.21351276130636e-01 ...
             +5.06899908932230e-01 +5.87715757240762e-01 +6.63044266930215e-01 ...
             +7.32182118740289e-01 +7.94483795967942e-01 +8.49367613732570e-01 ...
             +8.96321155766052e-01 +9.34906075937740e-01 +9.64762255587506e-01 ...
             +9.85611511545268e-01 +9.97263861849481e-01];         
w       =   [+7.01861000947021e-03 +1.62743947309057e-02 +2.53920653092621e-02 ...
             +3.42738629130214e-02 +4.28358980222263e-02 +5.09980592623761e-02 ...
             +5.86840934785358e-02 +6.58222227763622e-02 +7.23457941088488e-02 ...
             +7.81938957870697e-02 +8.33119242269464e-02 +8.76520930044040e-02 ...
             +9.11738786957641e-02 +9.38443990808043e-02 +9.56387200792743e-02 ...
             +9.65400885147275e-02 +9.65400885147285e-02 +9.56387200792748e-02 ...
             +9.38443990808043e-02 +9.11738786957645e-02 +8.76520930044033e-02 ...
             +8.33119242269462e-02 +7.81938957870706e-02 +7.23457941088489e-02 ...
             +6.58222227763618e-02 +5.86840934785352e-02 +5.09980592623760e-02 ...
             +4.28358980222264e-02 +3.42738629130212e-02 +2.53920653092626e-02 ...
             +1.62743947309056e-02 +7.01861000946989e-03];  
%%
hm    	=   (b-a)/2;
hp     	=   (b+a)/2;
f       =   func(hm*x+hp);
I       =   hm*(f*w');
end
%%