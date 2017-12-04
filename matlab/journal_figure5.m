clear all;
clc;

%% ====== parameters ======
disableIntraBoundaryFilter = 0;
width = 8;

subplot(1, 2, 1);

arrCG_resi_copy = [];
arrCG_resi_opt_error = [];
arrCG_resi_opt = [];
for sigma = 0 : 0.05 : 0.5
    N = width ^ 2;
    N_ext = N + width * 4 + 1;
    rho = 0.99; % correlation
    alpha = 0; % direction
    eta = 5; % strength of the directionality
    predModeIntra = 10;

    %% copying-based
    % ====== generate covariance matrix ======
    [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, disableIntraBoundaryFilter);
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
    
    
    %% optimal w/ error
    % ====== get optimal pred_mtx ======
    opt_pred_mtx_error = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights_error = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, i);
        opt_pred_mtx_error(i, ref_indices) = opt_weights_error;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx_error) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx_error)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt_error = [arrCG_resi_opt_error CG];
    
    
    %% optimal w/o error
    [pred_mtx, cov_mtx_ext_no_error] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, 0, disableIntraBoundaryFilter);
    % ====== get optimal pred_mtx ======
    opt_pred_mtx = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights = getOptimalPredictionWeights(cov_mtx_ext_no_error, ref_indices, i);
        opt_pred_mtx(i, ref_indices) = opt_weights;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
end

h0 =  plot(0 : 0.05 : 0.5, arrCG_resi_opt_error, '-vb'); hold on;
h1 = plot(0 : 0.05 : 0.5, arrCG_resi_opt, '-^k'); hold on;
h2 = plot(0 : 0.05 : 0.5, arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (26)', 'Optimal weights in Eq. (23)', 'Copying-based');
h_xl = xlabel('\sigma (reference deviation)'); 
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
grid on;


%% second figure

subplot(1, 2, 2);

arrCG_resi_copy = [];
arrCG_resi_opt_error = [];
arrCG_resi_opt = [];
for sigma = 0 : 0.05 : 0.5
    N = width ^ 2;
    N_ext = N + width * 4 + 1;
    rho = 0.99; % correlation
    alpha = -pi/4; % direction
    eta = 5; % strength of the directionality
    predModeIntra = 18;

    %% copying-based
    % ====== generate covariance matrix ======
    [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, disableIntraBoundaryFilter);
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_copy = [arrCG_resi_copy CG];
    
    
    %% optimal w/ error
    % ====== get optimal pred_mtx ======
    opt_pred_mtx_error = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights_error = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, i);
        opt_pred_mtx_error(i, ref_indices) = opt_weights_error;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx_error) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx_error)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt_error = [arrCG_resi_opt_error CG];
    
    
    %% optimal w/o error
    [pred_mtx, cov_mtx_ext_no_error] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, 0, disableIntraBoundaryFilter);
    % ====== get optimal pred_mtx ======
    opt_pred_mtx = zeros(size(pred_mtx));
    ref_indices = getHevcRefIdx( width, predModeIntra );
    for i = (N_ext - N + 1) : N_ext
        opt_weights = getOptimalPredictionWeights(cov_mtx_ext_no_error, ref_indices, i);
        opt_pred_mtx(i, ref_indices) = opt_weights;
    end
    
    % ====== apply prediction matrix ======
    resi_cov_mtx_ext = (eye(N_ext) - opt_pred_mtx) * cov_mtx_ext * (eye(N_ext) - opt_pred_mtx)';
    resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
    CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
    arrCG_resi_opt = [arrCG_resi_opt CG];
end

h0 =  plot(0 : 0.05 : 0.5, arrCG_resi_opt_error, '-vb'); hold on;
h1 = plot(0 : 0.05 : 0.5, arrCG_resi_opt, '-^k'); hold on;
h2 = plot(0 : 0.05 : 0.5, arrCG_resi_copy, 'r'); hold on;
h_legend = legend('Optimal weights in Eq. (26)', 'Optimal weights in Eq. (23)', 'Copying-based');
h_xl = xlabel('\sigma (reference deviation)'); 
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
grid on;

