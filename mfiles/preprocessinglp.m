function [P,Pi] = preprocessinglp(t,x,q,tau,ni,Ntheta,Nrho)
%normalize coordinates
[Nt,Nx]=size(t);[Nq,Ntau]=size(q);
T0=t(1);T1=t(end)+t(2)-t(1);X1=x(end)+x(1,2)-x(1,1);
t=(t)/(T1);x=x/X1;
q=q*X1/(T1);tau=(tau)/(T1);

%lp parameters
q2=q.^2;gamma=0;
alpha=(max(atan(q2(:)))+min(atan(q2(:))))/2;
beta=(max(atan(q2(:)))-min(atan(q2(:))));
aR=sqrt(2)*sin(beta)*sqrt(1/(2 + 2*sin(2*alpha)*(cos(beta) + sin(beta)) + sin(2*beta)));
Ox= -sqrt(0.2e1) * sin(alpha + pi / 0.4e1) * sin(beta) * aR / (-0.1e1 + cos(beta)) / 0.2e1;
Oy= sqrt(0.2e1) * aR * cos(alpha + pi / 0.4e1) * (0.1e1 - cos(beta)) / sin(beta) / 0.2e1;
recx=[-1/2 1/2 1/2 -1/2]*aR;
recy=[-1/2 -1/2 1/2 1/2]*aR;
[recx,recy]=mrotate(recx,recy,-alpha);
recx=recx+Ox;recy=recy+Oy;
y0=aR*cos(alpha-gamma)/cos(gamma)+recy(1);
x0=-aR*sin(alpha-gamma)/cos(gamma)+recx(1);
k=1/tan(beta/2);
b=y0-k*x0;
xm1=b/(-tan(beta/2)-k);
ym1=-tan(beta/2)*xm1;
am=sqrt(xm1^2+ym1^2);
am=am+aR*(T0/T1)^2*0.5;%decrease interval size due to the factor T0/T1

dtheta=(2*beta)/Ntheta;drho=-log(am)/Nrho;
thetas=linspace(-beta,beta,Ntheta+1);thetas=thetas(1:end-1);
rhos=linspace(log(am),0,Nrho+1);rhos=rhos(1:end-1);
[theta,rho]=ndgrid(thetas,rhos);

%t,x in lp
tt1=(t.^2-1/2)*aR;xx1=(x.^2-1/2)*aR;
[tt1,xx1]=mrotate(tt1,xx1,-alpha);
tt1=tt1+Ox;xx1=xx1+Oy;
thetat=atan2(xx1,tt1);rhot=log(sqrt(tt1.^2+xx1.^2));
emul=exp(reshape(rhot,size(t)));
thetat=thetat/2/beta;rhot=(rhot(:)-log(am)/2)/(-log(am));thetat=thetat(:);rhot=rhot(:);
%q,tau in lp
thetaq=-atan(q.^2)+alpha;
ttx=aR*(tau.^2-0.5)*cos(alpha)+aR*sin(alpha)*(1/2)+Ox;
tty=aR*(tau.^2-0.5)*sin(alpha)-aR*cos(alpha)*(1/2)+Oy;
rhoq=log((ttx+tty.*tan(thetaq)).*cos(thetaq));
cosmul=cos(thetaq-alpha);
thetaq=thetaq/2/beta;rhoq=(rhoq(:)-log(am)/2)/(-log(am));thetaq=thetaq(:);rhoq=rhoq(:);

%Jacobian hyperbolic integration, divided by 2x
J=4.*abs((aR.^2.*t)./(2.*(Ox.^2 + Oy.^2) + aR.^2.*(1 - 2.*t.^2 + 2.*t.^4 - 2.*x.^2 + 2.*x.^4) + 2.*aR.*(Ox.*(-1 + 2.*t.^2) + Oy.*(-1 + 2.*x.^2)).*cos(alpha) + 2.*aR.*(Ox - Oy + 2.*Oy.*t.^2 - 2.*Ox.*x.^2).*sin(alpha)));

%FT of zeta function

ker=@(x)((1/6*(2-abs(x)).^3).*(abs(x)<2).*(abs(x)>=1)+(2/3-1/2*abs(x).^2.*(2-abs(x))).*(abs(x)<1));%need [-N/2:N/2-1]
B3theta=ker(theta(:,1)/2/beta*Ntheta);B3rho=ker((rho(1,:)-log(am)/2)/log(am)*Nrho);B3thrho=B3theta*B3rho;
fB3=fftshift(fft2(ifftshift(B3thrho)));
fZ=fzeta_loop_weights(Ntheta,Nrho,Ntheta*dtheta,Nrho*drho,0,8);
fZB3=fftshift(fZ./fB3./fB3);
fZgpu=fZB3(1:end/2+1,:);fZgpu=[real(fZgpu(:)) imag(fZgpu(:))]';fZgpu=fZgpu(:);
fZintel=fZB3(1:end/2+1,:);fZintel=conj(fZintel);fZintel=[real(fZintel(:)');imag(fZintel(:)')];
%final normalization coefficient
cc=X1*(t(2,1)-t(1,1))*(x(1,2)-x(1,1))/(theta(2,1)-theta(1,1))/(rho(1,2)-rho(1,1))/aR;
J=J*cc;
%%sorting
ttn=thetat*Ntheta;xxn=rhot*Nrho;
add=2;
idthetatx=floor(ttn)+Ntheta/2+1+add;idrhotx=floor(xxn)+Nrho/2+1+add;
thetae=zeros(size(theta,1)+2*add,1);rhoe=zeros(size(rho,2)+2*add,1);thetae(add+1:end-add)=theta(:,1)/2/beta*Ntheta;rhoe(add+1:end-add)=(rho(1,:)-log(am)/2)/(-log(am))*Nrho;
dthetatx=ttn(:)-thetae(idthetatx);
drhotx=xxn(:)-rhoe(idrhotx);
idg=(idrhotx-1)*(Ntheta+2*add)+idthetatx;
[idg,reorids]=sort(idg);
st=ones(Ntheta+2*add,Nrho+2*add);nvals=zeros(Ntheta+2*add,Nrho+2*add);
k=1;
while(k<=Nt*Nx)
    j=1;
    st(idg(k))=k;    
    while((k+j<=Nt*Nx)&&(idg(k+j)==idg(k)))
        j=j+1;
    end;
    nvals(idg(k))=j;
    st(idg(k)+1:idg(min(k+j,Nt*Nx)))=k+j;
    k=k+j;
end
st(idg(Nt*Nx)+1:end)=st(idg(Nt*Nx))+j;%norm(st(idg+1)-st(idg)-nvals(idg))

%%fwd2
ttn=thetaq*Ntheta;
xxn=rhoq*Nrho;

idthetatauq=floor(ttn)+Ntheta/2+1+add;
idrhotauq=floor(xxn)+Nrho/2+1+add;

thetae=zeros(size(theta,1)+2*add,1);rhoe=zeros(size(rho,2)+2*add,1);thetae(add+1:end-add)=theta(:,1)/2/beta*Ntheta;rhoe(add+1:end-add)=(rho(1,:)-log(am)/2)/(-log(am))*Nrho;
dthetatauq=ttn(:)-thetae(idthetatauq);
drhotauq=xxn(:)-rhoe(idrhotauq);
%%%%%%
%%sorting adj
idgadj=(idrhotauq-1)*(Ntheta+2*add)+idthetatauq;
[idgadj,reoridsadj]=sort(idgadj);
stadj=ones(Ntheta+2*add,Nrho+2*add);nvalsadj=zeros(Ntheta+2*add,Nrho+2*add);
k=1;
while(k<=Nq*Ntau)
    j=1;
    stadj(idgadj(k))=k;    
    while((k+j<=Nq*Ntau)&&(idgadj(k+j)==idgadj(k)))
        j=j+1;
    end;
    nvalsadj(idgadj(k))=j;
    stadj(idgadj(k)+1:idgadj(min(k+j,Nq*Ntau)))=k+j;
    k=k+j;
end
stadj(idgadj(Nq*Ntau)+1:end)=stadj(idgadj(Nq*Ntau))+j;%norm(st(idg+1)-st(idg)-nvals(idg))

dthetatx=dthetatx(reorids);
drhotx=drhotx(reorids);
idthetatx=idthetatx(reorids);
idrhotx=idrhotx(reorids);

dthetatauq=dthetatauq(reoridsadj);
drhotauq=drhotauq(reoridsadj);
idthetatauq=idthetatauq(reoridsadj);
idrhotauq=idrhotauq(reoridsadj);

%to structure
Pi=[];%integer parameters
Pi.Nt=Nt;Pi.Nx=Nx;Pi.Nq=Nq;Pi.Ntau=Ntau;Pi.Ntheta=Ntheta;Pi.Nrho=Nrho;Pi.ni=ni;
Pi.st=st;Pi.nvals=nvals;Pi.reorids=reorids;
Pi.stadj=stadj;Pi.nvalsadj=nvalsadj;Pi.reoridsadj=reoridsadj;
Pi.idthetatx=idthetatx;Pi.idrhotx=idrhotx;
Pi.idthetatauq=idthetatauq;Pi.idrhotauq=idrhotauq;
P=[];%float parameters
P.t=t;P.x=x;P.q=q;P.tau=tau;P.theta=theta;P.rho=rho;
P.alpha=alpha;P.beta=beta;P.Ox=Ox;P.Oy=Oy;P.aR=aR;P.am=am;
P.J=J;P.fZ=fZ;P.fZB3=fZB3;P.fZintel=fZintel;P.fZgpu=fZgpu;P.cc=cc;P.fB3=fB3;
P.emul=emul;P.thetat=thetat;P.rhot=rhot;
P.cosmul=cosmul;P.thetaq=thetaq;P.rhoq=rhoq;
P.dthetatx=dthetatx;P.drhotx=drhotx;
P.dthetatauq=dthetatauq;P.drhotauq=drhotauq;


Pi=convert_toint(Pi);
P=convert_tofloat(P);
end

function [P] = convert_tofloat(P0)
P=[];
ff = fieldnames(P0);
    for i = 1:numel(ff)
    P.(ff{i})=single(P0.(ff{i}));  
    end
end

function [P] = convert_toint(P0)
P=[];
ff = fieldnames(P0);
    for i = 1:numel(ff)
    P.(ff{i})=int32(P0.(ff{i}));  
    end
end

function [ y1,y2 ] = mrotate( x1,x2,alpha )
t1=x1*cos(alpha)+x2*sin(alpha);
t2=-x1*sin(alpha)+x2*cos(alpha);
y1=t1(:);y2=t2(:);
end

function [dtheta,drho]=take_lpsteps(dt,dx,Ox,Oy,alpha,aR)
dtheta(1)=abs(atan((-2.*Oy + aR.*(-1 - 2.*(-2 + dx).*dx).*cos(alpha) + aR.*sin(alpha))./(-2.*Ox + aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dx).*dx).*sin(alpha))) + atan((2.*Oy + aR.*cos(alpha) - aR.*sin(alpha))./(-2.*Ox + aR.*(cos(alpha) + sin(alpha)))));
dtheta(2)=abs(atan((2.*Oy + aR.*cos(alpha) + aR.*(-1 + 2.*dt.^2).*sin(alpha))./(2.*Ox + aR.*(-1 + 2.*dt.^2).*cos(alpha) - aR.*sin(alpha))) + atan((2.*Oy + aR.*cos(alpha) - aR.*sin(alpha))./(-2.*Ox + aR.*(cos(alpha) + sin(alpha)))));
dtheta(3)=abs(atan((2.*Oy + aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dt).*dt).*sin(alpha))./(2.*Ox + aR.*(1 + 2.*(-2 + dt).*dt).*cos(alpha) - aR.*sin(alpha))) - atan((2.*Oy + aR.*(cos(alpha) + sin(alpha)))./(2.*Ox + aR.*cos(alpha) - aR.*sin(alpha))));
dtheta(4)=abs(atan((2.*Oy + aR.*(1 + 2.*(-2 + dx).*dx).*cos(alpha) + aR.*sin(alpha))./(-2.*Ox - aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dx).*dx).*sin(alpha))) + atan((2.*Oy + aR.*(cos(alpha) + sin(alpha)))./(2.*Ox + aR.*cos(alpha) - aR.*sin(alpha))));
dtheta(5)=abs(atan((Oy + (1./2).*aR.*((-1 + 2.*dx.^2).*cos(alpha) + sin(alpha)))./(Ox + (1./2).*aR.*(cos(alpha) + sin(alpha) - 2.*dx.^2.*sin(alpha)))) + atan(1 - (2.*(Ox + Oy + aR.*sin(alpha)))./(2.*Ox + aR.*(cos(alpha) + sin(alpha)))));
dtheta(6)=abs(atan((2.*Oy - aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dt).*dt).*sin(alpha))./(2.*Ox + aR.*(1 + 2.*(-2 + dt).*dt).*cos(alpha) + aR.*sin(alpha))) + atan(1 - (2.*(Ox + Oy + aR.*sin(alpha)))./(2.*Ox + aR.*(cos(alpha) + sin(alpha)))));
dtheta(7)=abs(atan((2.*Oy - aR.*(cos(alpha) + sin(alpha)))./(2.*Ox - aR.*cos(alpha) + aR.*sin(alpha))) + atan((-2.*Oy + aR.*(cos(alpha) + sin(alpha) - 2.*dt.^2.*sin(alpha)))./(2.*Ox + aR.*(-1 + 2.*dt.^2).*cos(alpha) + aR.*sin(alpha))));
dtheta(8)=abs(atan((2.*Oy - aR.*(cos(alpha) + sin(alpha)))./(2.*Ox - aR.*cos(alpha) + aR.*sin(alpha))) - atan((-2.*Oy + aR.*(cos(alpha) - 2.*dx.^2.*cos(alpha) + sin(alpha)))./(-2.*Ox + aR.*cos(alpha) + aR.*(-1 + 2.*dx.^2).*sin(alpha))));
dtheta=max(dtheta);

drho(1)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(-Ox + Oy).*cos(alpha) - aR.*(Ox + Oy).*sin(alpha)) - log((Oy + (1./2).*aR.*((1 + 2.*(-2 + dx).*dx).*cos(alpha) - sin(alpha))).^2 + (1./4).*(-2.*Ox + aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dx).*dx).*sin(alpha)).^2));
drho(2)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(-Ox + Oy).*cos(alpha) - aR.*(Ox + Oy).*sin(alpha)) - log((Ox + aR.*((-(1./2) + dt.^2).*cos(alpha) - sin(alpha)./2)).^2 + (1./4).*(2.*Oy + aR.*cos(alpha) + aR.*(-1 + 2.*dt.^2).*sin(alpha)).^2));
drho(3)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(Ox + Oy).*cos(alpha) + aR.*(-Ox + Oy).*sin(alpha)) - log((Ox + (1./2).*aR.*((1 + 2.*(-2 + dt).*dt).*cos(alpha) - sin(alpha))).^2 + (1./4).*(2.*Oy + aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dt).*dt).*sin(alpha)).^2));
drho(4)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(Ox + Oy).*cos(alpha) + aR.*(-Ox + Oy).*sin(alpha)) - log((Oy + (1./2).*aR.*((1 + 2.*(-2 + dx).*dx).*cos(alpha) + sin(alpha))).^2 + (Ox + (1./2).*aR.*(cos(alpha) + (-1 - 2.*(-2 + dx).*dx).*sin(alpha))).^2));
drho(5)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(Ox - Oy).*cos(alpha) + aR.*(Ox + Oy).*sin(alpha)) - log((Oy + (1./2).*aR.*((-1 + 2.*dx.^2).*cos(alpha) + sin(alpha))).^2 + (Ox + (1./2).*aR.*(cos(alpha) + sin(alpha) - 2.*dx.^2.*sin(alpha))).^2));
drho(6)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 + aR.*(Ox - Oy).*cos(alpha) + aR.*(Ox + Oy).*sin(alpha)) - log((1./4).*(2.*Oy - aR.*cos(alpha) + aR.*(1 + 2.*(-2 + dt).*dt).*sin(alpha)).^2 + (Ox + (1./2).*aR.*((1 + 2.*(-2 + dt).*dt).*cos(alpha) + sin(alpha))).^2));
drho(7)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 - aR.*(Ox + Oy).*cos(alpha) + aR.*(Ox - Oy).*sin(alpha)) - log((Ox + (1./2).*aR.*((-1 + 2.*dt.^2).*cos(alpha) + sin(alpha))).^2 + (1./4).*(-2.*Oy + aR.*(cos(alpha) + sin(alpha) - 2.*dt.^2.*sin(alpha))).^2));
drho(8)=(1./2).*abs(log(aR.^2./2 + Ox.^2 + Oy.^2 - aR.*(Ox + Oy).*cos(alpha) + aR.*(Ox - Oy).*sin(alpha)) - log((Oy + aR.*((-(1./2) + dx.^2).*cos(alpha) - sin(alpha)./2)).^2 + (Ox - (1./2).*aR.*(cos(alpha) + (-1 + 2.*dx.^2).*sin(alpha))).^2));
drho=max(drho);
end

