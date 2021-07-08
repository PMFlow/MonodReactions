%% Coupled degenerated flow & transport with Monod reactions
%% Deterministic-GRW flow solver; random-BGRW transport solver

clear all; close all
tic
global state;
initstate; state;
itest=1; % 0= loam soil; 1=clay soil;
%%   Grid Initialization
I=41; J=61; 
a=0; b=2;
c=0; d=3;
dx = (b-a)/(I-1);
x = a:dx:b;
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
y=c:dy:d;
y2=(y(1:J-1)+y(2:J))/2;
%%   Parameters
if itest ==1
    disp('Beit Netofa clay')
    Ksat = 8.2*10^-4;
    theta_res=0.0;
    theta_sat=0.446;
    alpha=0.152;
    nGM=1.17;
else
    disp('slit loam')
    Ksat = 4.96*10^-2; 
    theta_res=0.131;
    theta_sat=0.396;
    alpha=0.423;
    nGM=2.06;
end
T=1; % 3; % 5; % 
past=T/3;
t1=1; 
maxr=0.8;
D1=0.003; 
D2=D1;
Tolerance = 1e-5; 
S= 50000; 
Lp=1; Lc=1; 
%% Initialization of random hydraulic conductivity K
    Nmod = 100;
    ZC1 = 0.1;
    ZC2 = 0.01;
    varK = 4; 
    [X,Y] = meshgrid(x,y);
    [wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
    C1 = wavenum(:,1);
    C2 = wavenum(:,2);
    Ks= K_r(X(:)',Y(:)',Nmod,Ksat,varK,C1,C2,phi);
    Ks=reshape(Ks,J,I);
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p)  theta_GM(theta_res,theta_sat,p,alpha,nGM,ng);
K = @(tht) Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%% Initial Conditions
% IC - presure
p0 = 0.5 - 1*((y'-c)/(d-c))+0*x; % saturated/unsaturated regime ~ [Radu, 2004]
p = p0; pBC=p0;
% IC - concentrations
N=10^24;
c10 = zeros(J,I);
c20=5*ones(J,I);
i1=round(0.5/dx)+1; i2=round(1.5/dx)+1;
jb=round(b/6/dy);
i0=round(I/2); j0=round(0.9*J); 
dw=round(1/dx/10)+1;
c10(J,i1:i2)=2;
c20(J,i1:i2)=0;
c10=c10*N; c20=c20*N;
cBC1=c10; cBC2=c20;
c30=1e-3*ones(J,I);
c1=c10; c2=c20; c3=c30;
figure;
mesh(x,y,p);
xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
zlabel('$\psi(x,z,t=0)$','Interpreter','latex'); view(115,15);
figure;
mesh(x,y,c1/N);
xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
zlabel('$c_1(x,z,t=0)$','Interpreter','latex'); view(115,15);
figure;
mesh(x,y,c2/N);
xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
zlabel('$c_2(x,z,t=0)$','Interpreter','latex'); view(115,15);
%% Solution
tgraf=0; t=0; kt=1;
Vx=zeros(J,I); Vy=zeros(J,I); pp=zeros(J,I);
tht = theta(p); tht0=tht;
thtc10=tht0.*c10; thtc20=tht0.*c20;
pa=p; c1a=c10; c2a=c20; N1=N; N2=N1;
Dfactor=1.2;
dtc=Dfactor*(2*D1/Lc/dx^2+2*D2/Lc/dy^2); dtc=1/dtc;
convf=zeros(3,S); convc1=zeros(3,S); convc2=zeros(3,S);
while t<=T
    eps=zeros(1,S); epsc1=zeros(1,S); epsc2=zeros(1,S);
    DK=K(tht);
    Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
    Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
    D=DK(2:J-1,2:I-1);
    dtp=4*(max(max(Dx))/(Lp*dx^2)+max(max(Dy))/(Lp*dy^2)); dtp=maxr*1/dtp;
    dt=min(dtc,dtp);
    t=t+dt;
    %% Flow step
    for s=1:S
        DK=K(tht);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/Lp; ry=dt*Dy/dy^2/Lp;
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);
        %% Boundary conditions - flow
        %%%% BCY Left/Right
        pp(1:jb,1)=pBC(1:jb,1); % outflow on \Gmma_21
        pp(jb+1:J,1)=pp(jb+1:J,2); % no flux left
        pp(1:jb,I)=pBC(1:jb,I); % outflow on \Gmma_22 
        pp(jb+1:J,I)=pp(jb+1:J,I-1); % no flux right       
        %%%% BCX Bottom/Top
        pp(1,2:I-1)=pp(2,2:I-1)+dy; % no flux bottom
        if t<=t1 % linear p(t) on \Gamma_1
            pp(J,i1:i2)=-2.5+2.7*t/t1;
        else
            pp(J,i1:i2)=0.2; % constant \psi(t) on \Gamma_1
        end
        pp(J,2:i1-1)=pp(J-1,2:i1-1)-dy;  % no flux on top boundary - \Gamma_1
        pp(J,i2+1:I-1)=pp(J-1,i2+1:I-1)-dy;  % no flux on top boundary - \Gamma_1
        %% Source term - flow
        dtht=(tht0-tht)/Lp;
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1);
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        tht = theta(p);
        p0=p;
        %% Convergence criterion - flow
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        pa=p;
    end
    for s=1:S
        %% Transport step
        [Vx,Vy] = velocity(I,J,dx,dy,p0,D); % with gravity in Darcy's law
        thtc1=tht.*c1; thtc2=tht.*c2;
        nspec=2;
        for i=1:nspec
            switch i
                case 1
                    dthtc=(thtc10-thtc1)/Lc;
                    [c1t]=BGRW_2D_Monod_Richards_L(c10,cBC1,dthtc,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc,jb); % transport solver
                case 2
                    dthtc=(thtc20-thtc2)/Lc;
                    [c2t]=BGRW_2D_Monod_Richards_L(c20,cBC2,dthtc,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc,jb); % transport solver
            end
        end
        %% reaction
        [c1,c2,c3]=reaction_Monod3_L(c1t,c2t,c30,dt,tht,Lc,N1,N2);
        c10=c1; c20=c2; c30=c3;
        %% Convergence criterion - transport
        tol_epsc1=norm(c1/N1-c1a);
        tol_epsc2=norm(c2/N2-c2a);
        if kt*past>=t && kt*past<t+dt && t<=T
            epsc1(s)=tol_epsc1;
            epsc2(s)=tol_epsc2;
        end
        if max(tol_epsc1,tol_epsc2) <= Tolerance
            break
        end
        c1a=c1/N1; c2a=c2/N2;
    end
    if  kt*past>=t && kt*past<t+dt && t<=T
        tgraf=tgraf+1;
        rndt=kt*past;
        fprintf('t = %0.2e \n',rndt) ;
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        convf(kt,:)=eps;
        convc1(kt,:)=epsc1;
        convc2(kt,:)=epsc2;
        kt=kt+1;
    end
    tht0=tht; thtc10=tht.*c1; thtc20=tht.*c2;
    tE=t;
end
if  kt==4
    save('convf');
    save('convc1');
    save('convc2');
end
%% Results
kt_plot=kt-1;
plot_conv_Monod(kt_plot,S,strvect);
subplot_surf_c(x,y,c1/N,c2/N,c3);
subplot_surf_ptht(x,y,p,tht,Vx,Vy)
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)));
Pe=V*dx/D1;
maxV=max(max(sqrt(Vx.^2+Vy.^2)));
maxPe=maxV*dx/D1;
fprintf('meanV = %0.4e meanPe = %0.4e \n',V,Pe);
fprintf('maxV  = %0.4e maxPe  = %0.4e \n',maxV,maxPe);

toc %  (Laptop 16 GB, 2.11 GHz)

%% loam soil:
%% T=1
% The space step is : 5.00e-02 
% The time step is : 8.68e-05 
% meanV = 7.0544e-03 meanPe = 1.1757e-01 
% maxV  = 2.0526e-01 maxPe  = 3.4210e+00 
% Elapsed time is 3.3518 min.
% save('loamT1','x','y','c1','c2','c3','p','tht','Vx','Vy')
% iPe=sqrt(Vx.^2+Vy.^2)*dx/D1;
% [row,col] = find(iPe>2)
% row = 59 59 59
% col = 23 24 26
%% T=3
% The space step is : 5.00e-02 
% The time step is : 8.68e-05 
% meanV = 7.0415e-03 meanPe = 1.1736e-01 
% maxV  = 1.2979e-01 maxPe  = 2.1632e+00 
% Elapsed time is 7.0554 min.
% save('loamT3','x','y','c1','c2','c3','p','tht','Vx','Vy')
% row = 59
% col = 24
%% T=5
% The space step is : 5.00e-02 
% The time step is : 8.68e-05 
% meanV = 7.0639e-03 meanPe = 1.1773e-01 
% maxV  = 1.1271e-01 maxPe  = 1.8785e+00 
% Elapsed time is 10.0894 min.
% save('loamT5','x','y','c1','c2','c3','p','tht','Vx','Vy')

%% clay soil:
%% T=1:
% The space step is : 5.00e-02 
% The time step is : 5.25e-03 
% meanV = 7.9151e-05 meanPe = 1.3192e-03 
% maxV  = 3.7048e-03 maxPe  = 6.1746e-02 
% Elapsed time is 23.120326 seconds.
% save('clayT1','x','y','c1','c2','c3','p','tht','Vx','Vy')
%% T=3:
% The space step is : 5.00e-02 
% The time step is : 5.25e-03 
% meanV = 7.7958e-05 meanPe = 1.2993e-03 
% maxV  = 2.3710e-03 maxPe  = 3.9516e-02 
% Elapsed time is 57.160636 seconds.
% save('clayT3','x','y','c1','c2','c3','p','tht','Vx','Vy')
%% T=5:
% The space step is : 5.00e-02 
% The time step is : 5.25e-03 
% meanV = 7.6733e-05 meanPe = 1.2789e-03 
% maxV  = 2.2265e-03 maxPe  = 3.7108e-02 
% Elapsed time is 77.349145 seconds.
% save('clayT5','x','y','c1','c2','c3','p','tht','Vx','Vy')
