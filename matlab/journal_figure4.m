clear all;
clc;

%% ====== parameters ======
disableIntraBoundaryFilter = 1;

subplot(1, 2, 1);

arrCG_resi_copy = [];
arrCG_resi_opt_error = [];
arrCG_resi_opt = [];
for width = [4 8 16 32]
    
    N = width ^ 2;
    N_ext = N + width * 4 + 1;
    rho = 0.99; % correlation
    alpha = 0; % direction
    eta = 5; % strength of the directionality
    sigma = 0.0; % quantization error
    predModeIntra = 10;

    %% copying-based
    % ====== generate covariance matrix ======
    [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, disableIntraBoundaryFilter);
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
    
    
    %% optimal
    % ====== get optimal pred_mtx ======
    opt_pred_mtx = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, i);
        opt_pred_mtx(i, ref_indices) = opt_weights;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
end

h1 = plot([4 8 16 32], arrCG_resi_opt, '-^k'); hold on;
h2 = plot([4 8 16 32], arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (23)', 'Copying-based');
h_xl = xlabel('Block Size'); 
h_yl = ylabel('Coding gain (dB)');
h_t = title(['\rho=', sprintf('%4.2f', rho), ', \eta = 5, \alpha = 0 (mode 10)']);
grid on;

set(h1, 'LineWidth', 1.5);
set(h2, 'LineWidth', 1.5);
set(h_legend, 'FontSize', 15);
set(h_xl, 'FontSize', 15);
set(h_yl, 'FontSize', 15);
set(h_t, 'FontSize', 15);
set(gca,'FontSize', 12);
set(gca,'XTick',[4 8 16 32] );
grid on;


%% second figure

subplot(1, 2, 2);
arrCG_resi_copy = [];
arrCG_resi_opt_error = [];
arrCG_resi_opt = [];
for width = [4 8 16 32]
    
    N = width ^ 2;
    N_ext = N + width * 4 + 1;
    rho = 0.99; % correlation
    alpha = -pi/4; % direction
    eta = 5; % strength of the directionality
    sigma = 0.0; % quantization error
    predModeIntra = 18;

    %% copying-based
    % ====== generate covariance matrix ======
    [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, disableIntraBoundaryFilter);
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
    
    
    %% optimal
    % ====== get optimal pred_mtx ======
    opt_pred_mtx = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, i);
        opt_pred_mtx(i, ref_indices) = opt_weights;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
end

h1 = plot([4 8 16 32], arrCG_resi_opt, '-^k'); hold on;
h2 = plot([4 8 16 32], arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (23)', 'Copying-based');
h_xl = xlabel('Block Size'); 
h_yl = ylabel('Coding gain (dB)');
h_t = title(['\rho=', sprintf('%4.2f', rho), ', \eta = 5, \alpha = -\pi/4 (mode 18)']);
grid on;

set(h1, 'LineWidth', 1.5);
set(h2, 'LineWidth', 1.5);
set(h_legend, 'FontSize', 15);
set(h_xl, 'FontSize', 15);
set(h_yl, 'FontSize', 15);
set(h_t, 'FontSize', 15);
set(gca,'FontSize', 12);
set(gca,'XTick',[4 8 16 32] );
grid on;