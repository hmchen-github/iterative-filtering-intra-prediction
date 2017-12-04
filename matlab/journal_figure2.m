clear;
clc;

% ====== parameters ======
rho = 0.95;
N = 8;
sigma = 0.0; % quantization noise
arrCG_resi_opt = [];
arrCG_resi_opt_error = [];
arrCG_resi_copy = [];

% optimal
for sigma = 0 : 0.05 : 0.5
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
    
    % ======= prediction weight (optimal w/o error) =======
    pred_mtx( 2 : N + 1, 1 ) = -rho.^(1 : N);
    resi_cov_mtx_ext = pred_mtx * cov_mtx_ext * pred_mtx';
    resi_cov_mtx = resi_cov_mtx_ext( 2 : N + 1, 2 : N + 1 );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
    
    
    pred_mtx(1, 1) = 0;
    % ======= prediction weight (optimal w/ error) =======
    pred_mtx( 2 : N + 1, 1 ) = -rho.^(1 : N)/(1+sigma^2);
    
    resi_cov_mtx_ext = pred_mtx * cov_mtx_ext * pred_mtx';
    resi_cov_mtx = resi_cov_mtx_ext( 2 : N + 1, 2 : N + 1 );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt_error = [arrCG_resi_opt_error CG];
    
    % copying-based
    % ====== apply prediction matrix ======
    pred_mtx(1, 1) = 0;
    
    % ======= prediction weight (optimal) =======
    pred_mtx( 2 : N + 1, 1 ) = -1;
    
    resi_cov_mtx_ext = pred_mtx * cov_mtx_ext * pred_mtx';
    resi_cov_mtx = resi_cov_mtx_ext( 2 : N + 1, 2 : N + 1 );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
end

h0 = plot(0 : 0.05 : 0.5, arrCG_resi_opt_error, '-vb'); hold on;
h1 = plot(0 : 0.05 : 0.5, arrCG_resi_opt, '-^k'); hold on;
h2 = plot(0 : 0.05 : 0.5, arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (14)', 'Optimal weights in Eq. (6)', 'Copying-based');
h_xl = xlabel('\sigma (reference deviation)'); 
h_yl = ylabel('Coding gain (dB)');
grid on;

set(h0, 'LineWidth', 1.5);
set(h1, 'LineWidth', 1.5);
set(h2, 'LineWidth', 1.5);
set(h_legend, 'FontSize', 15);
set(h_xl, 'FontSize', 15);
set(h_yl, 'FontSize', 15);
set(gca,'FontSize', 12);
grid on;