function Px=compute_px(A,B,N,wl,tau)

K=size(A,1);
Nt=length(tau);
w_tilde=sqrt((A+wl).^2+B.^2);
alp=w_tilde*tau'; % K x Nt
bet=tau*wl; %Nt x 1
cos_alp=cos(alp);
sin_alp=sin(alp);
cos_bet=cos(bet);
sin_bet=sin(bet);

mz=B./w_tilde; % Kx1
mx=(A+wl)./w_tilde; %K x1

Phi=acos(cos_alp.*(ones(K,1)*cos_bet')- sin_alp.*(mx*sin_bet')); % K x Nt

M=1- (mz*ones(1,Nt)).^2.*((1-cos_alp).*(1-ones(K,1)*cos_bet')./(1+cos(Phi))).*sin(N*Phi/2).^2;

Px=(1-1e-6+prod(M,1)')/2;


