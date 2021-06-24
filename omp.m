% Returns a pair of lists of hyperfines A and B of size num_spins * 1.
% This pair of lists is ordered by highest correlation of the A,B pair's
% modeled signal with the signal to fit.
% A dictionary of spins sampled from the specified prior provides candidate
% spins.
% The prior can be 'uniform' (default) or 'grid'.
% This function will correlate signals in the specified domain: 'tau'
% (default), 'fourier_mag', 'db2'.

function [A, B] = omp(Y,M,N,wl,tau,num_spins,dict_size,prior,domain,residual_mode,a_min,a_max,b_min,b_max)    
    arguments
        Y (:,1) double
        M (1,1) double
        N (1,1) double
        wl (1,1) double
        tau (:,1) double
        num_spins (1,1) double
        dict_size (1,1) double
        prior (1,1) string = 'uniform'
        domain (1,1) string = 'tau'
        residual_mode (1,1) string = 'divide'
        a_min (1,1) double = -1e6
        a_max (1,1) double = 1e6
        b_min (1,1) double = -1e6
        b_max (1,1) double = 1e6
    end
    
    A = zeros(num_spins, 1);
    B = zeros(num_spins, 1);
    D = zeros(dict_size, 2);
    
    if strcmp(prior,'grid')
        % Create an orthogonal grid
        n = sqrt(dict_size);
        if n ~= floor(n)
            n = floor(n);
            dict_size = n^2;
            fprintf("%s %d %s\n", "Dictionary size changed to", dict_size, "to work with grid prior");
        end
        a = linspace(a_min, a_max, n);
        b = linspace(b_min, b_max, n);
        [A_grid, B_grid] = meshgrid(a,b);
        D = [reshape(A_grid,[size(A_grid,1)*size(A_grid,2),1]), reshape(B_grid,[size(B_grid,1)*size(B_grid,2),1])];
            
%             AA=ones(100,1)*linspace(-1e5,1e5,100);
%             D(:,1)=AA(:);
%             BB=(ones(100,1)*linspace(1e4,1e5,100))';
%             D(:,2)=BB(:);

    elseif strcmp(prior,'uniform')
        % Create dictionary D of random samples of the uniform prior
        D = rand(dict_size, 2);
        D(:,1) = a_min + 2*a_max * D(:,1);
        D(:,2) = b_min + 2*b_max * D(:,2);
%         D(:,2) = (2*(rand(dict_size,1)>0.5)-1).*(1e4+9.9e4*D(:,2));
%         D(:,2) = 2e3 + 1e4 * D(:,2);
    
    else
        fprintf("%s", "Specified prior not recognised");
    end
    
    if contains(domain, 'db')
        [cA, cD] = dwt(compute_px(D(1,1), D(1,2), N, wl, tau), domain);
        H = zeros(size(cA,1),dict_size);
        for k=1:dict_size
            [H(:,k), ~] = dwt(compute_mx(D(k,1), D(k,2), N, wl, tau), domain);
        end
        [e,~] = dwt(2*Y./M-1, domain);
        
        
    elseif strcmp(domain, 'fourier_mag')
        H = zeros(size(tau,1), dict_size);
        for k=1:dict_size
            % Magnitude of fourier transform of each spin's tau response
            H(:,k) = real(fft(compute_mx(D(k,1), D(k,2), N, wl, tau)));
        end
        % Initialise residual as magnitude of fourier transform of tau domain signal
        e = real(fft(2*Y./M-1));
        
    elseif strcmp(domain, 'tau')
        % One signal for every spin in the dictionary
        H = zeros(size(tau,1), dict_size);
        for k=1:dict_size
            H(:,k) = compute_mx(D(k,1), D(k,2), N, wl, tau);
        end
        % Initialise residual (tau domain signal)
        e = 2*Y./M-1;
        
    else
        fprintf("%s", "Specified domain not recognised");
    end

    Y_mean = mean(2*Y./M-1);
    % Total sum of squares
    SS_tot = sum( ( (2*Y./M-1)-Y_mean ).^2 );
    
    for k=1:num_spins
        % Find atom with maximum inner product
        [~, argmax] = min(sum((H-e*ones(1,size(H,2))).^2,1));
%         [~, argmax] = max(1-SS_tot./ ((e./H).^2));
        % TODO: Regularisation? euclidian distance between maximiser of
        % inner product and other spins?

        % Add corresponding (A,B) pair to sorted list of spins
        A(k) = D(argmax,1)
        B(k) = D(argmax,2)

        % Update residual by subtracting atom
        if strcmp(residual_mode,'divide')
            e = e ./ H(:,argmax);
        elseif strcmp(residual_mode,'subtract')
            e = e - H(:,argmax);
        end
        
        % Residual sum of squares
        SS_res = sum(e.^2);
        % R squared coefficient of determination
%         R_sq = 1 - SS_tot/SS_res
        R_sq = 1- exp(log(SS_tot)-log(SS_res))
        
        sum(e)

        % Remove atom from H
        H(:,argmax) = [];

        % Remove corresponding spin from dict
        D(argmax,:) = [];
    end
