function [A B MM]=death(Y,M,MM,A,B,N,wl,tau,lam,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 

%% MH for A
% Generate candidate using a random walk with variance sa
K=length(A);

i0=randi(K);
MM_cand=MM;
MM_cand(:,i0)=[];
P_cand=(1-1e-6+prod(MM_cand,2))/2;
P=(1-1e-6+prod(MM,2))/2;


% Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
plik=compute_LR(Y,P,P_cand,M,s2) -log(lam)+K;

% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A(i0)=[];
    B(i0)=[];
    MM=MM_cand;
end
