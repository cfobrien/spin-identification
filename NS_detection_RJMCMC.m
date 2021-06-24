function [ba bb A_map B_map bpost]=NS_detection_RJMCMC(A, B, K, Y,M0,N,wl,tau,Niter)

% K=1; %number of spins

% Initialization of the buffers to monitore the Markov chains
ba=zeros(1,Niter);
bb=zeros(1,Niter);

lam=15; % mean of Poisson distribution
% You can change the init. of the points to see the impact on the
% convergence

% A=5e3*ones(K,1);
% B=A;

% A=2*pi*[7.9e3] +1e4*randn; % We are cheating here and initialise close to the ground truth
% B=2*pi*[1e3]+1e4*randn;

M=M0;
sM=10;
% A=-1e6+2e6*rand(K,1);
% B=-1e6+2e6*rand(K,1);
% variances for the Gaussian random walks
sa=5e2;
sb=10*sa;
s2 = 20;% s2=1e0;
s0=1e7;
s_split=1e5;
A_map=A;
B_map=B;
Post_map=-inf;
for k=1:K
    MM(:,k)=compute_mx(A(k),B(k),N,wl,tau);
end
% TT=50; % Length of the sliding windows
% Q=rand(TT,1);
% Va=zeros(TT,1);
% Va(Q>0.55)=1;
% Vb=Va;

for m_compt=1:Niter
    if mod(m_compt,1000)==0
    [m_compt K M]
    figure(1)
    subplot(311)
%     plot(ba(1,1:m_compt-1)')
%     hold on
%     plot(bb(1,1:m_compt-1)')
%     hold off
plot(Y,'r')
hold on
plot(M*P,'b')
    hold off
    legend('Data','Fit')
    subplot(323)
    plot(bpost(1:m_compt-1))
    subplot(324)
    plot(log(bs2(1:m_compt-1)))
    subplot(325)
    plot(A,B,'+')
    subplot(326)
    plot(bK(1:m_compt-1))
    pause(0.01)
    end
    
    
    % Sequential update of the spin positions (MH)
    II=randi(K,min(10,K));
    for k=length(II)
        [A B MM conva convb]=update_AkBk(Y,M,MM,II(k),A,B,N,wl,tau,s2);
%         Va=[ Va(2:end); conva];
%         Vb=[ Vb(2:end); convb];
    end
    
    % Sequential update of the spin positions (MH)
    II=randi(K,min(10,K));
    for k=length(II)
        [A B MM]=update_AkBk_unif(Y,M,MM,k,A,B,N,wl,tau,s2);
    end
    

    
%     if m_compt>50
        
        u=rand;
        if u<0.05
            % Death
            if K>1
                [A B MM]=death(Y,M,MM,A,B,N,wl,tau,lam,s2);
            end
        elseif u>0.05 && u<0.1
            % Birth
            if K<70
            [A B MM]=birth(Y,M,MM,A,B,N,wl,tau,lam,s2);
            end
        elseif u>0.1 && u<0.55
            if K<70
            % Split
            [A B MM]=split(Y,M,MM,A,B,N,wl,tau,lam,s_split,s2);
            end
        else
            % Merge
            if K>1
            [A B MM]=merge(Y,M,MM,A,B,N,wl,tau,lam,s_split,s2);
            end
        end
%     end
    
    K=size(A,1);
    P=(1-1e-6+prod(MM,2))/2;
    if m_compt>1e4
    del2=1./(P'*P/s2+1/sM);
    mu=del2*(Y'*P/s2+M0/sM);
    M=mu+sqrt(del2)*randn;
    end
    s2=1./gamrnd(length(Y)/2,2./sum((Y-M*P).^2));
    
%     bpost(m_compt)=Y'*log(P)+(M-Y)'*log(1-P)+K*log(lam)-log(factorial(K))+2*K*log(1/2e5);
    bpost(m_compt)=-length(Y)/2*log(s2)-0.5/s2*sum((Y-M*P).^2)+K*log(lam)-log(factorial(K))+2*K*log(1/2e6);
    
    if bpost(m_compt)>Post_map
        A_map=A;
        B_map=B;
        Post_map=bpost(m_compt);
    end
    ba(1,m_compt)=A(1);
    bb(1,m_compt)=B(1);
    bs2(m_compt)=s2;
    bK(m_compt)=K;
end

