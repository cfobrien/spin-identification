function [A B MM]=split(Y,M,MM,A,B,N,wl,tau,lam,s_split,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 

%% MH for A
% Generate candidate using a random walk with variance sa
% a_cand=sqrt(s0)*randn;
% b_cand=sqrt(s0)*randn;
K=length(A);
i0=randi(K);
u=sqrt(s_split)*randn;
v=sqrt(s_split)+randn;
a_cand1=A(i0)+u;
a_cand2=A(i0)-u;
b_cand1=B(i0)+v;
b_cand2=B(i0)-v;
m_cand1=compute_mx(a_cand1,b_cand1,N,wl,tau)';
m_cand2=compute_mx(a_cand2,b_cand2,N,wl,tau)';
% if min(m_cand)>-1.1 && max(m_cand)<=1.1

MM_cand=MM;
MM_cand(:,i0)=m_cand1;
MM_cand=[MM m_cand2];
P_cand=(1-1e-6+prod(MM_cand,2))/2;
P=(1-1e-6+prod(MM,2))/2;


% Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P)) +2/2*log(2*pi*s_split)+u^2/(2*s_split)+v^2/(2*s_split)+log(4)+2*log(1/(2e5))+log(lam)-(K+1);
plik=compute_LR(Y,P,P_cand,M,s2)+2/2*log(2*pi*s_split)+u^2/(2*s_split)+v^2/(2*s_split)+log(4)+2*log(1/(2e6))+log(lam)-(K+1);
% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A=[A;a_cand2];
    B=[B;b_cand2];
    A(i0)=a_cand1;
    B(i0)=b_cand1;
    MM=MM_cand;
end

% end