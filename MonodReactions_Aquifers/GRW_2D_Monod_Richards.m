function [c]=GRW_2D_Monod_Richards(c0,cBC,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,stepU,dw)
%2D unbiased GRW algorithm for advection-diffusion processes;

%% GRW solution
vv=(Vy)*dt/dx./tht; uu=(Vx)*dt/dy./tht;
vr=vv; ur=uu; %V-approximated near the boundary
u=floor(ur*stepU+0.5);
v=floor((vr)*stepU+0.5);
% v=floor((vr+1)*stepU+0.5); % +1 ---> zero_velocity ! (because vr=-1 !!!)
d=2; d1=d; d2=d; % % -for both original (D=0.0001) and contours_numdiff(D=0.001, smaller dx)
rx=2*D1*dt/(d^2*dx^2)./tht; ry=2*D2*dt/(d^2*dy^2)./tht; r=rx+ry; %rx+ry<=1 !
n=c0; nBC=cBC;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; restf=0;
if max(max(r))>1
    r
    return
end

for y=(d1+stepU+1):(J-d1-stepU-1)
    for x=(d2+stepU+1):(I-d2-stepU-1)
        
        if n(y,x) > 0
            xa=x+u(y,x); ya=y+v(y,x);
            restr=n(y,x)*(1-r)+restr; nsta=floor(restr);
            restr=restr-nsta; njump=n(y,x)-nsta;
            if ya<1 || ya>J
                ya=y;
            end
            if xa<1 || xa>I
                xa=x;
            end
            if ya<1
                ya=1;
            end
            if ya>J
                ya=J;
            end
            if xa<1 
                xa=1;
            end
            if xa>I
                xa=I;
            end
            nn(ya,xa)=nn(ya,xa)+nsta;
            restjump=njump*rx/r+restjump;
            njumpx=floor(restjump); restjump=restjump-njumpx;
            njumpy=njump-njumpx;
            if(njumpy)>0
                restjumpy=njumpy/2+restjumpy;
                nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                nj(2)=njumpy-nj(1);
                if ya<d+1
                    nn(ya+d,xa)=nn(ya+d,xa)+nj(2); 
                elseif ya>J-d
                    nn(ya-d,xa)=nn(ya-d,xa)+nj(1); 
                else
                    for i=1:2
                        yd=ya+(2*i-3)*d;
                        nn(yd,xa)=nn(yd,xa)+nj(i);
                    end
                end
            end
            if(njumpx)>0
                restjumpx=njumpx/2+restjumpx;
                nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                nj(2)=njumpx-nj(1);
                if xa<d+1
                    nn(ya,xa+d)=nn(ya,xa+d)+nj(2); 
                elseif xa>I-d
                    nn(ya,xa-d)=nn(ya,xa-d)+nj(1); 
                else
                    for i=1:2
                        xd=xa+(2*i-3)*d;
                        nn(ya,xd)=nn(ya,xd)+nj(i);
                    end
                end
            end
        end
    end
end
%% boundary conditions
%%%% BCXBottom/Top
nn(:,1:j0-dw-1)=nBC(:,1:j0-dw-1); % set to IC on bottom boundary %% J=>x; I=>y !!!
nn(:,I-j0+dw+1:I)=nBC(:,I-j0+dw+1:I); % set to IC on top boundary
%%%% BCYLeft/Right
    nn(1:j0-dw-1,2:I-1)=nBC(1:j0-dw-1,2:I-1);  %  set to IC on left boundary
nn(j0-dw:j0+dw,i0-dw:i0+dw)=nBC(j0-dw:j0+dw,i0-dw:i0+dw);  % injection

n=nn; c=n;
