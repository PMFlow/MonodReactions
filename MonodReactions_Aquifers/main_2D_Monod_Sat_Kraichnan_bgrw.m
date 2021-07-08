%% Coupled 2D flow & reactive transport in saturated aquifers:
%% Kraichnan velocities; random BGRW transport solver; Monod reactions

clear all; close all
tic
global state;
initstate; state;
irand=1; % 0=constant Ksat; 1=random Ksat;
%% Grid Initialization
I=601; J=2401;
a=0; b=20;
c=0; d=80;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
%% Transport parameters
T= 30; % 300; 
Ksat = 25;
tht = 1;
D1=0.01; D2=D1;
stepU=2; 
p0 = 0.8 - 0.8*((y'-c)/(d-c))+0*x; % fully saturated regime ~ [Cirpka et al., 1999]
U_MEAN=-Ksat*((p0(J,2)-p0(1,2))/(d-c))/tht;
Dfactor=1.2;
dtc=Dfactor*(2*D1/dx^2/tht+2*D2/dy^2/tht); dt=1./dtc; 
%% Initialization Karaichnan routine
NMOD=100; % number of random modes
K_MEAN=Ksat; % meand K-field
varK=0.5; % variance of lnK-field
ZC1=2.0; ZC2=2.0; % correlation lengths
lambda=0.0; % width of spatial filtering
[wavenum, phi, amplitude] = V_Kraichnan_Gauss_param(NMOD,varK,ZC1,ZC2,U_MEAN,lambda);
%% Computation Kraichnan velocity
Vx=zeros(J,I); Vy=zeros(J,I);
for j=1:J
    for i=1:I
        if irand==0
            ur=0;
            vr=0;
        else
            yy=j*dy; xx=i*dx;
            [ur,vr] = V_Kraichnan_Gauss_func(xx,yy,wavenum, phi, amplitude);
        end
        Vx(j,i)=ur;
        Vy(j,i)=vr+U_MEAN;
    end
end
%% Initial Conditions - concentrations
N=10^24;
c10 = zeros(J,I);
c20 = 5*ones(J,I);
i0=round(I/2); j0=round(2.5/dx); dw=round(1/dx);
c10(j0-dw:j0+dw,i0-dw:i0+dw)=2; % injection in a square centerd at (j0,i0)
c20(j0-dw:j0+dw,i0-dw:i0+dw)=0;
c10=c10*N; c20=c20*N;  
cBC1=c10; cBC2=c20; 
c30=1e-3*ones(J,I);
c1=c10; c2=c20; c3=c30; 
figure;
mesh(x,y,c1/N);
xlabel('z','Interpreter','latex'); ylabel('x','Interpreter','latex');
zlabel('$c_1(x,z,t=0)$','Interpreter','latex'); view(115,15);
figure;
mesh(x,y,c2/N);
xlabel('z','Interpreter','latex'); ylabel('x','Interpreter','latex');
zlabel('$c_2(x,z,t=0)$','Interpreter','latex'); view(115,15);
%% Solution of the reactive transport
N1=N; N2=N1;
t=0; 
while t<=T
    t=t+dt;
    nspec=2; % mobile species
    for i=1:nspec
        switch i
            case 1
                [c1t]=BGRW_2D_Monod_Richards(c10,cBC1,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,dw); % transport solver
            case 2
                [c2t]=BGRW_2D_Monod_Richards(c20,cBC2,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,dw); % transport solver
        end
    end
    [c1,c2,c3]=reaction_Monod3(c1t,c2t,c30,dt,I,J);
    c10=c1; c20=c2; c30=c3;
end
%% Results
plot_fig_Monod_Sat(x,y,c1/N,c2/N,c3,Vx,Vy)
subplot_surf(x,y,c1/N,c2/N,c3);

fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)));
Pe=V*dx/D1;
maxV=max(max(sqrt(Vx.^2+Vy.^2)));
maxPe=maxV*dx/D1;
fprintf('meanV = %0.4e meanPe = %0.4e \n',V,Pe);
fprintf('maxV  = %0.4e maxPe  = %0.4e \n',maxV,maxPe);
toc % (Laptop 16 GB, 2.11 GHz)
%% variance=0.5; correlation lenths=2; T=300: meanPe=0.9365; maxPe=1.7507; (dx=0.0333; dt=0.0231; U=0.25)
% save('bgrwT300','x','y','c1','c2','c3') % Elapsed time is 1.2646 hours.
% save('bgrwT30','x','y','c1','c2','c3') % Elapsed time is 6.2937 min.