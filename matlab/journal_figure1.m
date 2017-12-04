clear;
clc;

% ====== parameters ======
rho = 0.95;
sigma = 0.0; % quantization noise
arrCG_resi_opt = [];
arrCG_resi_copy = [];

% optimal
for N = 4 : 32
    % ====== generate covariance matrix ======
    cov_mtx_ext = zeros( N + 1, N + 1 );
    cov_mtx = zeros( N, N );
    pred_mtx = eye( N + 1 );
    resi_cov_mtx_ext = zeros( N + 1, N + 1 );
    resi_cov_mtx = zeros( N, N );
    for i = 1 : N + 1
        for j = 1 : N + 1
            cov_mtx_ext(i, j) = rho^(abs(i - j));
        end
    end
    cov_mtx_ext(1, 1) = cov_mtx_ext(1, 1) + sigma ^ 2;
    
    % ====== apply prediction matrix ======
    pred_mtx(1, 1) = 0;
    
    % ======= optimal prediction weight =======
    pred_mtx( 2 : N + 1, 1 ) = -rho.^(1 : N)/(1+sigma^2);
    
    resi_cov_mtx_ext = pred_mtx * cov_mtx_ext * pred_mtx';
    resi_cov_mtx = resi_cov_mtx_ext( 2 : N + 1, 2 : N + 1 );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
    
    % copying-based
    % ====== apply prediction matrix ======
    pred_mtx(1, 1) = 0;
    
    % ======= optimal prediction weight =======
    pred_mtx( 2 : N + 1, 1 ) = -1;
    
    resi_cov_mtx_ext = pred_mtx * cov_mtx_ext * pred_mtx';
    resi_cov_mtx = resi_cov_mtx_ext( 2 : N + 1, 2 : N + 1 );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
end


h1 = plot(4 : 32, arrCG_resi_opt, '-^k'); hold on;
h2 = plot(4 : 32, arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (6)', 'Copying-based');
h_xl = xlabel('N (Number of samples)'); 
h_yl = ylabel('Coding gain (dB)');
grid on;

set(h1, 'LineWidth', 1.5);
set(h2, 'LineWidth', 1.5);
set(h_legend, 'FontSize', 15);
set(h_xl, 'FontSize', 15);
set(h_yl, 'FontSize', 15);
set(gca,'FontSize', 12);
grid on;