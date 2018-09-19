addpath('mfiles','mex_gpu');
N=512;
Nx=N;Nt=N;Nq=N;Ntau=N;T0=0;T1=4;X1=5;
t=T0+(0:Nt-1)/Nt*(T1-T0);x=(0:Nx-1)/Nx*X1;
q=0.3+(0:Nq-1)/Nq*0.6;tau=T0+(0:Ntau-1)/Ntau*(T1-T0);
[t,x]=ndgrid(t,x);[q,tau]=ndgrid(q,tau);

%hyperbolas
f=zeros(size(t));
rickh=@(t,cent,a)(1-2*pi^2*a^2*(t).^2).*exp(-pi^2*a^2*t.^2); 
ql=[0.6 0.6 0.6 0.6 0.6 0.5 0.45 0.7]*max(q(:));tl=[0.1 0.3 0.5 0.65 0.8 0.3 0.5 0.5]*(T1-T0)+T0;
for l=1:8, for j=1:size(x,2),f(:,j)=f(:,j)+rickh(t(:,1)-sqrt(tl(l).^2+ql(l)^2*x(1,j).^2),tl(l),6);end;end;f=f*(t(2,1)-t(1,1));
%smooth near t=T.
v=ones(Nt,Nx);v(end-15:end,:)=repmat(cos(linspace(0,1,16)*pi/2)',1,Nx);
f=f.*v;

%size of the log-polar grid 
Ntheta=N;Nrho=N;
[P,Pi] = preprocessinglp(t,x,q,tau,1,Ntheta,Nrho);
cl_lp=create_class_gpu(Pi,P);

f=single(f);
%forward transform
Rlp=fwd(cl_lp,f);
%adjoint transform
flp=adj(cl_lp,Rlp);

%adjoint test
(sum(flp(:).*f(:))-sum(Rlp(:).*Rlp(:)))/sum(Rlp(:).*Rlp(:))
delete(cl_lp)