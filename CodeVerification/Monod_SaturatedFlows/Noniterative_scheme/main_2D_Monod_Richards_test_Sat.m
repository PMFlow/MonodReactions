%% Deterministic-GRW flow solver; Noniterative random-BGRW solver for transport and Monod reactions

clear all; close all
tic

Itest = 1;
L2_p = zeros(4,1); eoc_p=zeros(3,1);
L2_c1 = zeros(4,1); eoc_c1=zeros(3,1);
L2_c2 = zeros(4,1); eoc_c2=zeros(3,1);
 for it = [10 20]% 40  80] 
%%   Grid Initialization
    I = it+1; J = 3*it/2+1;
    a=0; b=2;
    c=0; d=3;
    dx = (b-a)/(I-1);
    x = a:dx:b;
    dy=(d-c)/(J-1);
    y=c:dy:d;
    [X,Y] = meshgrid(x,y) ;    
%%  Parameters
    T=1;
    past=T/3;
    S= 50000;
    maxr=1; 
    Tolerance = 1e-6;
    Lp=1;
    D1=0.025; D2=D1;
%%  Exact solutions and sources
    solEp = @(t,x,y)  t.*x.*(2-x).*y.*(3-y); % exact presure solution for saturated regime
    solEc1 = @(t,x,y) x.*(2-x).*y.^3*exp(-0.1*t)/27;
    solEc2 = @(t,x,y) (x-1).^2.*y.^2*exp(-0.1*t)/9;    
    K = @(p) 5e-2*ones(size(p));
%%   Initial Conditions 
% IC - pressure  
    p0 = solEp(0,X,Y);
    p = p0; pBC=p0;
% IC concentration
    N=10^24; N1=N; N2=N1;
    c10=solEc1(0,X,Y)*N; c1=c10;
    c20=solEc2(0,X,Y)*N; c2=c20;
%% Solution
tgraf=0; t=0; kt=1; 
Vx=zeros(J,I); Vy=zeros(J,I);
convf=zeros(3,S); pp=zeros(J,I);
tht = theta(p);
tht0=tht; 
pa=p; 
Dfactor=2;
dtc=Dfactor*(2*D1/dx^2/min(min(tht))+2*D2/dy^2/min(min(tht))); dtc=1/dtc; % BGRW
while t<=T
    DK=K(p);
    Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
    Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
    D=DK(2:J-1,2:I-1);
    [fs,fg,fh,VxE,VyE]=F_saturated(t,x,y);
    dtc=Dfactor*(2*D1/dx^2/min(min(tht))+2*D2/dy^2/min(min(tht))); dtc=1/dtc; 
    dtp=4*(max(max(Dx))/(Lp*dx^2)+max(max(Dy))/(Lp*dy^2)); dtp=maxr*1/dtp;
    dt=min(dtc,dtp);
    t=t+dt;
    eps=zeros(1,S);     
%% Flow step        
    for s=1:S
        DK=K(p);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);        
        rx=dt*Dx/dx^2/Lp; ry=dt*Dy/dy^2/Lp;
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);       
        %% Boundary conditions - pressure
        %%%% BCYLeft/Right
        pp(:,1)=pBC(:,1); % Dirichlet, analytical solution
        pp(:,I)=pBC(:,I);
        %%%% BCXBottom/Top      
        pp(1,2:I-1)=pBC(1,2:I-1);
        pp(J,2:I-1)=pBC(J,2:I-1);
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1) ...
            + fs(2:J-1,2:I-1)*dt/Lp;
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        p0=p; 
        %% Convergence criterion (L-scheme for Richards eq.)
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        pa=p;     
    end
    if  kt*past>=t && kt*past<t+dt && t<=T
        tgraf=tgraf+1;
        rndt=kt*past;
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        convf(kt,:)=eps;
        kt=kt+1;
    end
%% Transport step
[Vx,Vy] = velocity(I,J,dx,dy,p0,D,DK); % with gravity in Darcy's law
        nspec=2;
        for i=1:nspec
            switch i
                case 1
                    cBC1=solEc1(t,X,Y)*N1;
                    dthtc1=(tht0-tht).*c1./tht0;
                    [c1t]=BGRW_2D_Monod_Richards(c1,cBC1,tht,I,J,dx,dy,dt,Vx,Vy,D1,D2); % transport solver
                case 2
                    cBC2=solEc2(t,X,Y)*N2;
                    dthtc2=(tht0-tht).*c2./tht0;
                    [c2t]=BGRW_2D_Monod_Richards(c2,cBC2,tht,I,J,dx,dy,dt,Vx,Vy,D1,D2); % transport solver
            end
        end
%% reaction
    [c1,c2]=reaction_Monod2(c1t,c2t,dt,I,J,N1,N2); % m~1e-3...
    c1=c1+dthtc1+fg./tht*N1*dt;
    c2=c2+dthtc2+fh./tht*N2*dt;
%    
    tht0=tht; 
    tE=t;
end
if  kt==4
        save('convf');
end

%% Results
    kt_plot=kt-1;
    plot_conv(kt_plot,S,strvect)
    fprintf('The space step is : %0.4e \n',dx) ;
    fprintf('The time step is : %0.4e \n',dt) ;
    fprintf('The total time is : %0.4e \n',tE) ;
    cE_p=solEp(tE,X,Y);
    L2_p(Itest) = ( dx * dy )^(1/2) *norm(p-cE_p); 
    c1E=solEc1(tE,X,Y); c1=c1/N1;
    L2_c1(Itest) = ( dx * dy )^(1/2) *norm(c1-c1E);
    c2E=solEc2(tE,X,Y); c2=c2/N2;
    L2_c2(Itest) = ( dx * dy )^(1/2) *norm(c2-c2E);
    if Itest >1
        eoc_p(Itest-1)=log10(L2_p(Itest-1)/L2_p(Itest))/log10(2);
        eoc_c1(Itest-1)=log10(L2_c1(Itest-1)/L2_c1(Itest))/log10(2);
        eoc_c2(Itest-1)=log10(L2_c2(Itest-1)/L2_c2(Itest))/log10(2);
    end
    Itest = Itest+1; 
 end
fprintf('L2_p  : %0.2e \n',L2_p)
fprintf('EOC_p : %0.2e \n',eoc_p)
fprintf('L2_c1  : %0.4e \n',L2_c1)
fprintf('EOC_c1 : %0.4e \n',eoc_c1)
fprintf('L2_c2  : %0.4e \n',L2_c2)
fprintf('EOC_c2 : %0.4e \n',eoc_c2)
%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)));
Pe=V*dx/D1;
maxV=max(max(sqrt(Vx.^2+Vy.^2)));
maxPe=maxV*dx/D1;
fprintf('meanV = %0.4e meanPe = %0.4e \n',V,Pe);
fprintf('maxV  = %0.4e maxPe  = %0.4e \n',maxV,maxPe);
toc %  (Laptop 16 GB, 2.11 GHz)

%% main_2D_Monod_Richards_test
% The space step is : 2.0000e-01 
% The time step is : 6.0000e-02 
% The total time is : 1.0200e+00 
% The space step is : 1.0000e-01 
% The time step is : 1.5000e-02 
% The total time is : 1.0050e+00 
% The space step is : 5.0000e-02 
% The time step is : 3.7500e-03 
% The total time is : 1.0013e+00 
% The space step is : 2.5000e-02 
% The time step is : 9.3750e-04 
% The total time is : 1.0003e+00 

% L2_p  : 1.76e-01 
% L2_p  : 4.44e-02 
% L2_p  : 1.21e-02 
% L2_p  : 7.22e-03 
                        % EOC_p : 1.99e+00 
                        % EOC_p : 1.87e+00 
                        % EOC_p : 7.50e-01 
% L2_c1  : 1.9174e-02 
% L2_c1  : 4.8337e-03 
% L2_c1  : 1.2667e-03 
% L2_c1  : 5.5664e-04 
                        % EOC_c1 : 1.9880e+00 
                        % EOC_c1 : 1.9321e+00 
                        % EOC_c1 : 1.1862e+00 
% L2_c2  : 2.5455e-02 
% L2_c2  : 5.8822e-03 
% L2_c2  : 1.4682e-03 
% L2_c2  : 4.8838e-04 
                        % EOC_c2 : 2.1135e+00 
                        % EOC_c2 : 2.0024e+00 
                        % EOC_c2 : 1.5879e+00 
% meanV = 1.1479e-01 meanPe = 2.2959e-01 
% maxV  = 2.2440e-01 maxPe  = 4.4879e-01 
% Elapsed time is 23.942375 seconds (for it=[10 20]).
% Elapsed time is 4.8297 hours  (for it=[10 20 40 80]).

% s between 800 and 4200 iterations of the GRW-flow L-scheme
