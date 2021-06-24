function [A B MM conva convb]=update_AkBk(Y,M,MM,k,A,B,N,wl,tau,s2)
% This function was originally designed for A and B vectors but in this
% demo it is used with A and B scalars.
% This codes performs 2 sequential updates for A and B using 1D MH updates 
conva=0;
convb=0;
sa=5e2;
sb=1*sa;
%% MH for A
% Generate candidate using a random walk with variance sa
a_cand=A(k)+sqrt(sa)*randn;
b_cand=B(k)+sqrt(sb)*randn;
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
% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P));

% If a prior is used, it should be included here. 
% Note that because the proposal is symmetric the ratio of the proposals
% cancels out.

% Accept/reject step 
if rand<exp(plik) % if accept
    A(k)=a_cand;
    B(k)=b_cand;
    P=P_cand;
    MM=MM_cand;
    conva=1;
end

% %% MH for A
% % Generate candidate using a random walk with variance sa
% a_cand=A(k)+sqrt(sa)*randn;
% 
% MM_cand=MM;
% m_cand=compute_mx(a_cand,B(k),N,wl,tau);
% MM_cand(:,k)=m_cand;
% P_cand=(1-1e-6+prod(MM_cand,2))/2;
% P=(1-1e-6+prod(MM,2))/2;
% % Compute likelihood for current point and candidate
% % P=compute_px(A,B,N,wl,tau);
% % P_cand=compute_px(A_cand,B,N,wl,tau);
% % Note there is no prior here: it is equavalent to using a uniform prior
% % distribution
% 
% % Compute log-ratio of the likelihoods (numerically more stable than computing the ratio directly)
% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P));
% % If a prior is used, it should be included here. 
% % Note that because the proposal is symmetric the ratio of the proposals
% % cancels out.
% 
% % Accept/reject step 
% if rand<exp(plik) % if accept
%     A(k)=a_cand;
%     P=P_cand;
%     MM=MM_cand;
%     conva=1;
% end
% 
% %% Same for B 
% b_cand=B(k)+sqrt(sb)*randn;
% 
% MM_cand=MM;
% m_cand=compute_mx(A(k),b_cand,N,wl,tau);
% MM_cand(:,k)=m_cand;
% P_cand=(1-1e-6+prod(MM_cand,2))/2;
% % P=(1-1e-6+prod(MM,1)')/2;
% 
% 
% % P_cand=compute_px(A,B_cand,N,wl,tau);
% 
% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P));
% 
% if rand<exp(plik)
%     B(k)=b_cand;
%     convb=1;
%     MM=MM_cand;
% end
% 
% 
% 
