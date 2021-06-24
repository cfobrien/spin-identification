function [A B MM]=merge(Y,M,MM,A,B,N,wl,tau,lam,s_split,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 

%% MH for A
% Generate candidate using a random walk with variance sa
K=length(A);
ord=randperm(K);
i0=ord(1);
i1=ord(2);
a_cand=(A(i0)+A(i1))/2;
u=(A(i0)-A(i1))/2;
b_cand=(B(i0)+B(i1))/2;
v=(B(i0)-B(i1))/2;
MM_cand=MM;
MM_cand(:,i0)=compute_mx(a_cand,b_cand,N,wl,tau)';
MM_cand(:,i1)=[];
P_cand=(1-1e-6+prod(MM_cand,2))/2;
P=(1-1e-6+prod(MM,2))/2;


% Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P))-2/2*log(2*pi*s_split)-u^2/(2*s_split)-v^2/(2*s_split)-log(4)-2*log(1/(2e5)) -log(lam)+K;
plik=compute_LR(Y,P,P_cand,M,s2)-2/2*log(2*pi*s_split)-u^2/(2*s_split)-v^2/(2*s_split)-log(4)-2*log(1/(2e6)) -log(lam)+K;
% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A(i0)=a_cand;
    B(i0)=b_cand;
    A(i1)=[];
    B(i1)=[];
    MM=MM_cand;
end
