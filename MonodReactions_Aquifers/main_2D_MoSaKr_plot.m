%% GRW solutions for different variances and correlation lengths of the velocity field
%% Kraichnan velocities; random unbiased GRW transport solver; Monod reactions

clear all; close all
global state;
tic
%%   Grid Initialization
I=161; J=641; 
a=0; b=20;
c=0; d=80;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
%%   Transport parameters
Ksat = 25;
tht=1;
T=300;
D=0.01; D1=D; D2=D;
N=10^24;
lambda=0.0; % width of spatial filtering
NMOD=100; % number of random modes
K_MEAN=Ksat; % meand K-field
stepU=2;
p0 = 0.8 - 0.8*((y'-c)/(d-c))+0*x;
U_MEAN=-Ksat*((p0(J,2)-p0(1,2))/(d-c))/tht;
dt=stepU*dx/abs(U_MEAN);
%% Initialization Karaichnan routine
figure;
for iZ=1:2
    if iZ==1
        kp=1;
    else
        kp=2;
    end
    for iK=1:3
        initstate; state;        
        if iZ==1
            ZC1=1.0; ZC2=1.0;
        else
            ZC1=2.0; ZC2=2.0;
        end
        if iK==1
            varK=0.1;
        elseif iK==2
            varK=0.5;
        else
            varK=1;
        end
        [wavenum, phi, amplitude] = V_Kraichnan_Gauss_param(NMOD,varK,ZC1,ZC2,U_MEAN,lambda);
        Vx=zeros(J,I); Vy=zeros(J,I);
        for j=1:J
            for i=1:I
                yy=j*dy; xx=i*dx;
                [ur,vr] = V_Kraichnan_Gauss_func(xx,yy,wavenum, phi, amplitude);
                Vx(j,i)=ur;
                Vy(j,i)=vr+U_MEAN;
            end
        end
        
        %%  Initial Conditions - concentrations 
        c10 = zeros(J,I);
        c20 = 5*ones(J,I); %
        i0=round(I/2); j0=round(2.5/dx); dw=round(1/dx);
        c10(j0-dw:j0+dw,i0-dw:i0+dw)=2; % injection in a square centerd at (j0,i0)
        c20(j0-dw:j0+dw,i0-dw:i0+dw)=0;
        c10=c10*N; c20=c20*N;
        cBC1=c10; cBC2=c20;
        c30=1e-3*ones(J,I);
        c1=c10; c2=c20; c3=c30;
        %% Solution
        t=0;
        while t<=T
            t=t+dt;
            nspec=2;
            for i=1:nspec
                switch i
                    case 1
                        [c1t]=GRW_2D_Monod_Richards(c10,cBC1,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,stepU,dw);
                    case 2
                        [c2t]=GRW_2D_Monod_Richards(c20,cBC2,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,stepU,dw);
                end
            end
            [c1,c2,c3]=reaction_Monod3(c1t,c2t,c30,dt,I,J);
            c10=c1; c20=c2; c30=c3;
        end
        
        %% Results
        c1=c1/N; c2=c2/N;
        L=num2str(ZC1); S=num2str(varK);
        
        subplot(9,2,kp)
        surf(y,x,c1','EdgeColor','none'); colorbar; view(0,90); caxis([0 2]);
        set(gca, 'XTick', []); set(gca, 'YTick', []); xlim([0 80]); 
        title(['$c_1(\lambda=$',L,'$,\sigma^2=$',S,')'],'Interpreter','latex','FontSize',7);
        kp=kp+2;
        subplot(9,2,kp)
        surf(y,x,c2','EdgeColor','none'); colorbar; view(0,90); caxis([0 5]);
        set(gca, 'XTick', []); set(gca, 'YTick', []); xlim([0 80]); 
        title(['$c_2(\lambda=$',L,'$,\sigma^2=$',S,')'],'Interpreter','latex','FontSize',7);
        kp=kp+2;
        subplot(9,2,kp)
        surf(y,x,c3','EdgeColor','none'); colorbar; view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []); xlim([0 80]); 
        title(['$c_3(\lambda=$',L,'$,\sigma^2=$',S,')'],'Interpreter','latex','FontSize',7);
        kp=kp+2;
    end
end
toc % Elapsed time is 23.059187 seconds (Laptop 16 GB, 2.11 GHz)
