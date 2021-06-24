%%
clear;
close all;
clc;

wl=2*pi*10.7084*1e6*40.3*1e-3; % 2 pi * parameter for carbon * magnetic field applied
N=256;
% M=5e3; %number of bernoulli trials per tau
M = 1;

tau = readNPY('tau_256_f.npy');

y = readNPY('y_256_f.npy');


ind = find(tau<25);
tau = tau(ind);
y = y(ind);

tau=tau*1e-6;

%% DNN
dnn_spins = csvread('dnn_spins.csv') .* (2 * pi * 1e3);
signal_dnn = compute_px(dnn_spins(:,1), dnn_spins(:,2), N, wl, tau);

dnn_snr = compare("DNN", signal_dnn, y, tau);

%% OMP
[A_omp, B_omp] = omp(y,M,N,wl,tau,30,10000,"grid","tau","divide");
% signal_omp = csvread('omp_spins_.csv');
signal_omp = compute_px(A_omp, B_omp, N, wl, tau);
omp_snr = compare("OMP", signal_omp, y, tau);

save results_omp.mat A_omp B_omp


%% RJMCMC (no OMP)
t=1;
% for t=1:10
%     t
Niter=5e6; %Number of MCMC iterations

K = 10;
A=-1e6+2e6*rand(K,1);
B=-1e6+2e6*rand(K,1);
[ba bb A_est B_est bpost]=NS_detection_RJMCMC(A, B, K, y*M,M,N,wl,tau,Niter);
AA{t}=A_est;
BB{t}=B_est;
Post(:,t)=bpost;
% end

save results_rjmcmc.mat AA BB Post
%%
load results_rjmcmc.mat
A_rjmcmc = AA{1};
B_rjmcmc = BB{1};
signal_rjmcmc = compute_px(A_rjmcmc, B_rjmcmc, N, wl, tau);
rjmcmc_snr = compare("RJMCMC", signal_rjmcmc, y, tau);


%% RJMCMC w/ OMP
t=1;
% for t=1:10
%     t
Niter=1e6; %Number of MCMC iterations

load results_omp.mat
K = size(A_omp,1);
A=A_omp;
B=B_omp;
[ba bb A_est B_est bpost]=NS_detection_RJMCMC(A, B, K, y*M,M,N,wl,tau,Niter);
AA{t}=A_est;
BB{t}=B_est;
Post(:,t)=bpost;
% end

save results_rjmcmc_omp.mat AA BB Post
%%
load results_rjmcmc_omp.mat
A_rjmcmc_omp = AA{1};
B_rjmcmc_omp = BB{1};
signal_rjmcmc_omp = compute_px(A_rjmcmc_omp, B_rjmcmc_omp, N, wl, tau);
rjmcmc_omp_snr = compare("RJMCMC w/ OMP", signal_rjmcmc_omp, y, tau);


%% RJMCMC w/ OMP + GAUS BLUR
t=1;
% for t=1:10
%     t
Niter=5e5; %Number of MCMC iterations

load results_omp.mat
K = size(A_omp,1);
%%
y_blur = gaussian_blur(y, 10, 5);

[ba bb A_est B_est bpost]=NS_detection_RJMCMC(A_omp, B_omp, K, y_blur*M,M,N,wl,tau,Niter);
AA{t}=A_est;
BB{t}=B_est;
Post(:,t)=bpost;
% end

save results_rjmcmc_omp_gaus_blur.mat AA BB Post
%%
load results_rjmcmc_omp_gaus_blur.mat
A_rjmcmc_omp_gaus_blur = AA{1};
B_rjmcmc_omp_gaus_blur = BB{1};
signal_rjmcmc_omp_gaus_blur = compute_px(A_rjmcmc_omp_gaus_blur, B_rjmcmc_omp_gaus_blur, N, wl, tau);
rjmcmc_omp_gaus_blur_snr = compare("RJMCMC w/ OMP + gaus blur", signal_rjmcmc_omp_gaus_blur, y, tau);


%% OMP plus blur preliminary run + rjmcmc
t=1;
% for t=1:10
%     t
Niter=5e5; %Number of MCMC iterations

load results_rjmcmc_omp_gaus_blur.mat
K = size(A_rjmcmc_omp_gaus_blur,1);
[ba bb A_est B_est bpost]=NS_detection_RJMCMC(A_rjmcmc_omp_gaus_blur, B_rjmcmc_omp_gaus_blur, K, y_blur*M,M,N,wl,tau,Niter);
AA{t}=A_est;
BB{t}=B_est;
Post(:,t)=bpost;
% end
save results_rjmcmc_blur_prelim.mat AA BB Post
%%
load results_rjmcmc_blur_prelim.mat
A_rjmcmc_blur_prelim = AA{1};
B_rjmcmc_blur_prelim = BB{1};
signal_rjmcmc_blur_prelim = compute_px(A_rjmcmc_blur_prelim, B_rjmcmc_blur_prelim, N, wl, tau);
rjmcmc_blur_prelim_snr = compare("RJMCMC w/ prelim blur", signal_rjmcmc_blur_prelim, y, tau);


%%
figure
xlabel("B (KHz)");
ylabel("A (KHz)");
hold ON
scatter(dnn_spins(:,2), dnn_spins(:,1));
scatter(B_omp, A_omp);
scatter(B_rjmcmc_omp, A_rjmcmc_omp);
scatter(B_rjmcmc, A_rjmcmc);
scatter(B_rjmcmc_omp_gaus_blur, A_rjmcmc_omp_gaus_blur);
legend("DNN", "OMP", "RJMCMC", "RJ+OMP", "RJ+OMP+BLUR");


%%
snr = @(signal_gt, signal_est) 10*log10((rms(signal_gt)/rms(signal_est - signal_gt))^2);

function s2n = compare(name, signal_est, y, tau)
    s2n = snr(y, signal_est);
    figure
    plot(tau, y);
    hold ON
    plot(tau, signal_est);
    legend(["Y" name]);
    title(strcat("SNR ", num2str(s2n), " dB"));
end


