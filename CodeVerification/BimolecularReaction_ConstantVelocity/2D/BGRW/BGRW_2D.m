function [c]=BGRW_2D(c0,c0E,fg,I,J,dx,dy,dt,Vx,Vy,D1,D2)
% BGRW function - transport step in coupled flow & trnasport 2D

%% BGRW solution
v=Vy*dt/dx; u=Vx*dt/dy;
ru=2*D1*dt/dx^2*ones(J,I);
rv=2*D2*dt/dy^2*ones(J,I);
n=c0; nBC=c0E;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; restf=0;
    for y=1:J
        for x=1:I 
            if n(y,x) > 0
                rx=ru(y,x); ry=rv(y,x); r=rx+ry;
                restr=n(y,x)*(1-r)+restr; nsta=floor(restr);
                restr=restr-nsta; njump=n(y,x)-nsta;
                nn(y,x)=nn(y,x)+nsta;
                restjump=njump*ry/r+restjump;
                njumpy=floor(restjump); restjump=restjump-njumpy;
                njumpx=njump-njumpy;
                if(njumpy)>0
                    restjumpy=njumpy*0.5*(1-v(y,x)/ry)+restjumpy;
                    nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                    nj(2)=njumpy-nj(1);
                    if y==1
                        nn(2,x)=nn(2,x)+nj(2); 
                        nn(1,x)=nn(1,x)+nj(1);
                    elseif y==J
                        nn(J-1,x)=nn(J-1,x)+nj(1); 
                        nn(J,x)=nn(J,x)+nj(2);
                    else
                        for i=1:2
                            yd=y+(2*i-3);
                            nn(yd,x)=nn(yd,x)+nj(i);
                        end
                    end
                end                
                if(njumpx)>0
                    restjumpx=njumpx*0.5*(1-u(y,x)/rx)+restjumpx;
                    nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                    nj(2)=njumpx-nj(1);
                    if x==1
                        nn(y,2)=nn(y,2)+nj(2); 
                        nn(y,1)=nn(y,1)+nj(1);
                    elseif x==I
                        nn(y,I-1)=nn(y,I-1)+nj(1); 
                        nn(y,I)=nn(y,I)+nj(2);
                    else
                        for i=1:2
                            xd=x+(2*i-3);
                            nn(y,xd)=nn(y,xd)+nj(i);
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
    n=nn;
    c=n;