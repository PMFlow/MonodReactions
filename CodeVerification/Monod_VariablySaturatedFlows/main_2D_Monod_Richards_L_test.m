%% Deterministic-GRW flow solver; Iterative random-BGRW L-scheme for transport and Monod reactions
%% test for degenerated Richards equation and transition from unsaturated to saturated flow regime

clear all; close all
tic

Itest = 1;
L2_p = zeros(4,1); eoc_p=zeros(3,1);
L2_Vx = zeros(4,1); eoc_Vx=zeros(3,1);
L2_Vy = zeros(4,1); eoc_Vy=zeros(3,1);
L2_c1 = zeros(4,1); eoc_c1=zeros(3,1);
L2_c2 = zeros(4,1); eoc_c2=zeros(3,1);
for it = [10 20]% 40 80] 
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
    Lp=1; Lc=1; 
    Dfactor=2; 
    D1=0.025; D2=D1;
    %%  Exact solutions and sources
    solEp = @(t,x,y)  -t.*x.*(2-x).*y.*(3-y)+x/4+y/4; % exact solution with transition from unsaturated to saturated flow regime
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
    pp=zeros(J,I);
    tht = theta(p); tht0=tht;  
    thtc10=tht0.*c10; thtc20=tht0.*c20;
    pa=p; c1a=c10/N; c2a=c20/N;
    dtc=Dfactor*(2*D1/Lc/dx^2+2*D2/Lc/dy^2); dtc=1/dtc; 
    convf=zeros(3,S); convc1=zeros(3,S); convc2=zeros(3,S);
    while t<=T
        eps=zeros(1,S); epsc1=zeros(1,S); epsc2=zeros(1,S);
        DK=K(p);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        [fs,fg,fh,VxE,VyE]=F(t,x,y);
        dtp=4*(max(max(Dx))/(Lp*dx^2)+max(max(Dy))/(Lp*dy^2)); dtp=maxr*1/dtp;
        dt=min(dtc,dtp);
        t=t+dt;
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
            %% Boundary conditions - flow
            %% Diriclet
            %%%% BCYLeft/Right
            pp(:,1)=pBC(:,1); % Dirichlet, analytical solution
            pp(:,I)=pBC(:,I); 
            %%%% BCXBottom/Top
            pp(1,2:I-1)=pBC(1,2:I-1); 
            pp(J,2:I-1)=pBC(J,2:I-1); 
            %% Source term - flow
            dtht=(tht0-tht)/Lp;
            f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1) ...
                + fs(2:J-1,2:I-1)*dt/Lp;
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
            [Vx,Vy] = velocity(I,J,dx,dy,p0,D,DK); % with gravity in Darcy's law
            thtc1=tht.*c1; thtc2=tht.*c2;
            nspec=2;
            for i=1:nspec
                switch i
                    case 1
                        cBC1=solEc1(t,X,Y)*N1;
                        dthtc1=(thtc10-thtc1)/Lc;
                        [c1t]=BGRW_2D_Monod_Richards_L(c1,cBC1,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc); % transport solver
                    case 2
                        cBC2=solEc2(t,X,Y)*N2;
                        dthtc2=(thtc20-thtc2)/Lc;
                        [c2t]=BGRW_2D_Monod_Richards_L(c2,cBC2,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc); % transport solver
                end
            end
            %% reaction
            [c1,c2]=reaction_Monod2_L(c1t,c2t,dt,tht,N1,N2,Lc);

            c1=c1+dthtc1+fg/Lc*N1*dt;
            c2=c2+dthtc2+fh/Lc*N2*dt;
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
    fprintf('The space step is : %0.4e \n',dx) ;
    fprintf('The time step is : %0.4e \n',dt) ;
    fprintf('The total time is : %0.4e \n',tE) ;
    cE_p=solEp(tE,X,Y);
    L2_p(Itest) = ( dx * dy )^(1/2) *norm(p-cE_p);
    L2_Vx(Itest) = ( dx * dy )^(1/2) *norm(Vx-VxE);
    L2_Vy(Itest) = ( dx * dy )^(1/2) *norm(Vy-VyE);
    c1E=solEc1(tE,X,Y); c1=c1/N1;
    L2_c1(Itest) = ( dx * dy )^(1/2) *norm(c1-c1E);
    c2E=solEc2(tE,X,Y); c2=c2/N2;
    L2_c2(Itest) = ( dx * dy )^(1/2) *norm(c2-c2E);
    if Itest >1
        eoc_p(Itest-1)=log10(L2_p(Itest-1)/L2_p(Itest))/log10(2);
        eoc_Vx(Itest-1)=log10(L2_Vx(Itest-1)/L2_Vx(Itest))/log10(2);
        eoc_Vy(Itest-1)=log10(L2_Vy(Itest-1)/L2_Vy(Itest))/log10(2);
        eoc_c1(Itest-1)=log10(L2_c1(Itest-1)/L2_c1(Itest))/log10(2);
        eoc_c2(Itest-1)=log10(L2_c2(Itest-1)/L2_c2(Itest))/log10(2);
    end
    Itest = Itest+1;
end
fprintf('L2_p  : %0.2e \n',L2_p)
fprintf('EOC_p : %0.2e \n',eoc_p)
fprintf('L2_Vx  : %0.4e \n',L2_Vx)
fprintf('EOC_Vx : %0.4e \n',eoc_Vx)
fprintf('L2_Vy  : %0.4e \n',L2_Vy)
fprintf('EOC_Vy : %0.4e \n',eoc_Vy)
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

%% main_2D_Monod_Richards_L_test.m
% The space step is : 2.0000e-01 
% The time step is : 1.0000e-01 
% The total time is : 1.1000e+00 
% The space step is : 1.0000e-01 
% The time step is : 2.5000e-02 
% The total time is : 1.0000e+00 
% The space step is : 5.0000e-02 
% The time step is : 6.2500e-03 
% The total time is : 1.0062e+00 
% The space step is : 2.5000e-02 
% The time step is : 1.5625e-03 
% The total time is : 1.0000e+00 
% L2_p  : 2.86e-01 
% L2_p  : 7.18e-02 
% L2_p  : 1.83e-02 
% L2_p  : 4.94e-03 
                                    % EOC_p : 2.00e+00 
                                    % EOC_p : 1.98e+00 
                                    % EOC_p : 1.88e+00 
% L2_Vx  : 1.6900e-02 
% L2_Vx  : 6.0169e-03 
% L2_Vx  : 2.2330e-03 
% L2_Vx  : 8.0407e-04 
                                    % EOC_Vx : 1.4900e+00 
                                    % EOC_Vx : 1.4300e+00 
                                    % EOC_Vx : 1.4736e+00 
% L2_Vy  : 6.0560e-03 
% L2_Vy  : 2.1824e-03 
% L2_Vy  : 8.0345e-04 
% L2_Vy  : 3.0958e-04 
                                    % EOC_Vy : 1.4725e+00 
                                    % EOC_Vy : 1.4416e+00 
                                    % EOC_Vy : 1.3759e+00 
                                    
% L2_c1  : 2.2940e-02 
% L2_c1  : 5.8149e-03 
% L2_c1  : 1.6672e-03 
% L2_c1  : 6.8327e-04 
                                    % EOC_c1 : 1.9801e+00 
                                    % EOC_c1 : 1.8023e+00 
                                    % EOC_c1 : 1.2869e+00 
% L2_c2  : 3.5801e-02 
% L2_c2  : 9.3713e-03 
% L2_c2  : 2.6333e-03 
% L2_c2  : 8.0962e-04 
                                    % EOC_c2 : 1.9337e+00 
                                    % EOC_c2 : 1.8314e+00 
                                    % EOC_c2 : 1.7016e+00 
% meanV = 1.1957e-01 meanPe = 1.1957e-01 
% maxV  = 2.4253e-01 maxPe  = 2.4253e-01 
% Elapsed time is 6.325224 seconds (for it=[10 20]).
% Elapsed time is 1.1676/1.1699 hours (for it=[10 20 40 80]).
