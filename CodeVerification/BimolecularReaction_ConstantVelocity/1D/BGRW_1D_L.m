function [c]=BGRW_1D_L(c0,c0E,dc,fg,I,dx,dt,Lc,q,D)
% BGRW function - transport step in coupled flow & trnasport 1D

%% BGRW solution
u=q*dt/Lc/dx;
ru=2*D*dt/Lc/dx^2*ones(1,I);
n=c0; n0E=c0E;
nn=zeros(1,I);
restr=0; restjump=0; restf=0;

for x=1:I
    if n(x) > 0
        r=ru(x);
        restr=n(x)*(1-r)+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        nn(x)=nn(x)+nsta;
        if(njump)>0
            restjump=njump*0.5*(1-u(x)/r)+restjump;
            nj(1)=floor(restjump); restjump=njump-nj(1);
            nj(2)=floor(restjump); restjump=restjump-nj(2);
            if x==1
                nn(2)=nn(2)+nj(2);
            elseif x==I
                nn(I-1)=nn(I-1)+nj(1);
            else
                for i=1:2
                    xd=x+(2*i-3);
                    nn(xd)=nn(xd)+nj(i);
                end
            end
        end
    end
end
%% boundary conditions
%%%% BC_Left/Right
nn(1)=n0E(1);
nn(I)=n0E(I);
%% Source term
restf=dc(2:I-1)/Lc+fg(2:I-1)*dt/Lc+restf; nf=floor(restf); restf=restf-nf;

nn(2:I-1)=nn(2:I-1)+nf;
n=nn;
c=n;