function [c]=GRW_2D(c0,cBC,fg,tht,I,J,dx,dy,dt,Vx,Vy,D1,D2,stepU)
%2D unbiased GRW algorithm for advection-diffusion processes

%% BGRW solution
vv=(Vy)*dt/dx/tht; uu=(Vx)*dt/dy/tht;
vr=vv; ur=uu; 
u=floor(ur*stepU+0.5);
v=floor((vr+1)*stepU+0.5); 
d=1;
rx=2*D1*dt/(d^2*dx^2)/tht; ry=2*D2*dt/(d^2*dy^2)/tht; r=rx+ry; %rx+ry<=1 
n=c0; nBC=cBC;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; restf=0;
if r>1
    r
    return
end

for y=1:J
    for x=1:I 
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
    %% Boundary conditions - concentration
    %%%% BCYLeft/Right
    nn(:,1)=nBC(:,1); 
    nn(:,I)=nBC(:,I);
    %%%% BCXBottom/Top
    nn(1,:)=nBC(1,:);
    nn(J,:)=nBC(J,:);
    %% Source term
    restf=fg(2:J-1,2:I-1)*dt+restf; nf=floor(restf); restf=restf-nf;
    nn(2:J-1,2:I-1)=nn(2:J-1,2:I-1)+nf;

n=nn; c=n;
c(1:J-1,:)=c(2:J,:);

