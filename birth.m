function [A B MM]=birth(Y,M,MM,A,B,N,wl,tau,lam,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 

%% MH for A
% Generate candidate using a random walk with variance sa
% a_cand=sqrt(s0)*randn;
% b_cand=sqrt(s0)*randn;
% a_cand=-1e5+2e5*rand;
% b_cand=-1e5+2e5*rand;
a_cand=-1e6+2e6*rand;
b_cand=-1e6+2e6*rand;
K=length(A);
m_cand=compute_mx(a_cand,b_cand,N,wl,tau)';
if min(m_cand)>-1.1 && max(m_cand)<=1.1
    
MM_cand=[MM m_cand];
P_cand=(1-1e-6+prod(MM_cand,2))/2;
P=(1-1e-6+prod(MM,2))/2;


% Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
plik=compute_LR(Y,P,P_cand,M,s2) +log(lam)-(K+1);
% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A=[A;a_cand];
    B=[B;b_cand];
    MM=MM_cand;
end

end