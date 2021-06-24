% clear
% close all
% clc

N=256;

wl=2*pi*10.7084*1e6*40.3*1e-3;

tau = readNPY('tau_256_f.npy');
y = readNPY('y_256_f.npy');

ind=find(tau<25); % all measurements between 10 and 10us tau
tau=tau(ind);
y=y(ind);
figure
plot(tau,y)

tau = tau ./ 1e6;

dnn_hyperfines = csvread('dnn_spins.csv') .* (2 * pi * 1e3);
px = compute_px(dnn_hyperfines(:,1), dnn_hyperfines(:,2), N, wl, tau);

figure
scatter(dnn_hyperfines(:,2), dnn_hyperfines(:,1));
xlabel("B (KHz)");
ylabel("A (KHz)");

figure
subplot(2,1,1);
plot(tau, y);
hold ON
plot(tau, px);
hold OFF
xlabel("tau (us)");
ylabel("Px");
legend("Y", "DNN");
subplot(2,1,2);
plot(tau(1:200), y(1:200));
hold ON
plot(tau(1:200), px(1:200));
xlabel("tau (us)");
ylabel("Px");
legend("Y", "DNN");