function snr=compute_snr(signal_gt, signal_est, N, wl, tau)
    noise = signal_est - signal_gt;
    snr = 10*log10((rms(signal_gt)/rms(noise))^2);