function plik=compute_LR(Y,P,P_cand,M,s2)


% plik=Y'*log(P_cand./P)+(M-Y)'*log((1-P_cand)./(1-P));
plik=-0.5/s2* ( sum((Y-M*P_cand).^2)-sum((Y-M*P).^2));