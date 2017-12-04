clear;
clc;

subplot(1, 2, 1);
% ====== parameters ======
width = 4;
N = width ^ 2;
N_ext = N + width * 4 + 1;
rho = 0.99; % correlation
alpha = 0; % direction
eta = 1; % strength of the directionality
sigma = 0.2;
arrIntraMode = 1;
iteration = 0 : 50;        
kernel = [0 1/6 0; 1/6 1/3 1/6; 0 1/6 0];

arrCG_resi = [];
arrCG_dct = [];
arrCG_dst = [];
arrCG_klt = [];

arrCG_filter_resi = [];

arrCG_opt_resi = [];

for iter = iteration
    
    cov_mtx_ext = zeros( N_ext, N_ext );
    cov_mtx = zeros( N, N );
    pred_mtx = zeros( N_ext, N_ext);
    resi_cov_mtx_ext = zeros( N_ext, N_ext );
    resi_cov_mtx = zeros( N, N );
    
    
    % ====== no filtering ======
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 0 );
        
        % ====== apply prediction matrix ======
        resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
        resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
        CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
        arrCG_resi = [arrCG_resi CG];
    end
    
    
    % ====== filtering ======
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 0 );

        filter_pred_mtx = filterPrediction( width, pred_mtx, iter, kernel );
        
        % ====== apply prediction matrix ======
        resi_cov_mtx_ext = (eye(N_ext) - filter_pred_mtx) * cov_mtx_ext * (eye(N_ext) - filter_pred_mtx)';
        resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
        CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
        arrCG_filter_resi = [arrCG_filter_resi CG];
    end
    
    % ====== optimal prediction ======    
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 1 );
        
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
        arrCG_opt_resi = [arrCG_opt_resi CG];
    end
end

plot(iteration, arrCG_opt_resi, 'k-^', 'LineWidth', 1.5); hold on;
plot(iteration, arrCG_filter_resi, 'k-.', 'LineWidth', 1.5); hold on;
plot(iteration, arrCG_resi, 'k', 'LineWidth', 1.5); hold on;

h_legend = legend('Optimal', 'Filtering', 'HEVC (copying-based)');
h_xlabel = xlabel('Iteration Number T');
h_ylabel = ylabel('Coding Gain (dB)');
h_title = title({'Coding gain comparison', [' ( 4\times 4, \rho = 0.99, \eta = ' num2str(eta) ', \alpha = 0, \sigma = ' num2str(sigma) ')']});

set(h_title, 'FontName','Times New Roman', 'FontSize', 15);
set(h_xlabel, 'FontName','Times New Roman', 'FontSize', 15);
set(h_ylabel, 'FontName','Times New Roman', 'FontSize', 15);
set(h_legend, 'FontName','Times New Roman', 'FontSize', 15);
grid on;


%% 8x8
subplot(1, 2, 2);
% ====== parameters ======
width = 8;
N = width ^ 2;
N_ext = N + width * 4 + 1;
rho = 0.99; % correlation
alpha = 0; % direction
eta = 1; % strength of the directionality
sigma = 0.2;
arrIntraMode = 1;
iteration = 0 : 5: 50;

arrCG_resi = [];
arrCG_dct = [];
arrCG_dst = [];
arrCG_klt = [];

arrCG_filter_resi = [];

arrCG_opt_resi = [];

for iter = iteration
    
    cov_mtx_ext = zeros( N_ext, N_ext );
    cov_mtx = zeros( N, N );
    pred_mtx = zeros( N_ext, N_ext);
    resi_cov_mtx_ext = zeros( N_ext, N_ext );
    resi_cov_mtx = zeros( N, N );
    
    
    % ====== no filtering ======
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 0 );
        
        % ====== apply prediction matrix ======
        resi_cov_mtx_ext = (eye(N_ext) - pred_mtx) * cov_mtx_ext * (eye(N_ext) - pred_mtx)';
        resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
        CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
        arrCG_resi = [arrCG_resi CG];
    end
    
    
    % ====== filtering ======
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 0 );

        filter_pred_mtx = filterPrediction( width, pred_mtx, iter, kernel );
        
        % ====== apply prediction matrix ======
        resi_cov_mtx_ext = (eye(N_ext) - filter_pred_mtx) * cov_mtx_ext * (eye(N_ext) - filter_pred_mtx)';
        resi_cov_mtx = resi_cov_mtx_ext( width * 4 + 2 : N_ext, width * 4 + 2 : N_ext );
        CG = -10*(sum(log10(diag(resi_cov_mtx))))/N;
        arrCG_filter_resi = [arrCG_filter_resi CG];
    end
    
    % ====== optimal prediction ======    
    for predModeIntra = arrIntraMode
        
        % ====== generate covariance matrix ======
        [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, 1 );
        
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
        arrCG_opt_resi = [arrCG_opt_resi CG];
    end
end

plot(iteration, arrCG_opt_resi, 'k-^', 'LineWidth', 1.5); hold on;
plot(iteration, arrCG_filter_resi, 'k-.', 'LineWidth', 1.5); hold on;
plot(iteration, arrCG_resi, 'k', 'LineWidth', 1.5); hold on;


h_legend = legend('Optimal', 'Filtering', 'HEVC (copying-based)');
h_xlabel = xlabel('Iteration Number T');
h_ylabel = ylabel('Coding Gain (dB)');
h_title = title({'Coding gain comparison', [' ( 8\times 8, \rho = 0.99, \eta = ' num2str(eta) ', \alpha = 0, \sigma = ' num2str(sigma) ')']});

set(h_title, 'FontName','Times New Roman', 'FontSize', 15);
set(h_xlabel, 'FontName','Times New Roman', 'FontSize', 15);
set(h_ylabel, 'FontName','Times New Roman', 'FontSize', 15);
set(h_legend, 'FontName','Times New Roman', 'FontSize', 15);
grid on;