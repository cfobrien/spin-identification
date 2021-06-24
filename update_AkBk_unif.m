function [A B MM]=update_AkBk_unif(Y,M,MM,k,A,B,N,wl,tau,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 


%% MH for A
% Generate candidate using a random walk with variance sa
% a_cand=-1e5+2e5*rand;
% b_cand=-1e5+2e5*rand;
a_cand=-1e6+2e6*rand;
b_cand=-1e6+2e6*rand;
MM_cand=MM;
m_cand=compute_mx(a_cand,b_cand,N,wl,tau);
MM_cand(:,k)=m_cand;
P_cand=(1-1e-6+prod(MM_cand,2))/2;
P=(1-1e-6+prod(MM,2))/2;
% Compute likelihood for current point and candidate
% P=compute_px(A,B,N,wl,tau);
% P_cand=compute_px(A_cand,B,N,wl,tau);
% Note there is no prior here: it is equavalent to using a uniform prior
% distribution

% Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
plik=compute_LR(Y,P,P_cand,M,s2);

% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A(k)=a_cand;
    B(k)=b_cand;
    P=P_cand;
    MM=MM_cand;
end