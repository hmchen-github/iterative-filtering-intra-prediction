clear;
clc;

% ====== parameters ======
rho = 0.99; % correlation along dominating direction
rho2 = rho^5; % correlation perpendicular to dominating direction
sigmaSet = [0 0.2 0.5]; % quantization error
etaList = [1 log(rho2)/log(rho)] % stength
isHEVC = 0;
modeSet = 10;
width = 4;
kernel = [0 1/6 0; 1/6 1/3 1/6; 0 1/6 0];

lineStyleSet = {'-k', '--r', ':b', '-*'};

for column = 1 : 4
    for row = 1 : 4
        posInVector = width * 4 + 1 + (column - 1) * width + row;
        
        subplot(4,4, row * 4 - 4 + column)
        
        weightsRecord = zeros(width * 4 + 1, size(sigmaSet, 1));
        count = 0;
        legendInfo = {};
        
        for sigma = sigmaSet
            count = count + 1;
            
            N = width ^ 2;
            N_ext = N + width * 4 + 1;
            
            % ====== optimal prediction ======
            for predModeIntra = modeSet
                
                % ====== mode-dependent parameters ======
                if (predModeIntra == 1)
                    eta = 1;
                    alpha = 0;
                else
                    eta = etaList(2);
                    alpha = (10 - predModeIntra) / 32 * pi;
                end

                % ====== generate covariance matrix ======
                [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( width, predModeIntra, rho, alpha, eta, sigma, ~isHEVC);
                
                % ====== filtering prediction ======
                filter_pred_mtx = filterPrediction( width, pred_mtx, 5, kernel );
                
                % permutate the reference pixels
                permutMtx = eye(N_ext);
                
                for i = 1 : 2 * width + 1
                    permutMtx(i, i) = 0;
                    permutMtx(i, 2 * width + 2 - i) = 1;
                end
                
                cov_mtx_ext = permutMtx * cov_mtx_ext * permutMtx';
                
                filter_pred_mtx = filter_pred_mtx * inv(permutMtx);
                
                % ====== get optimal pred_mtx (both no-error and error cases) ======
                opt_pred_mtx = zeros(size(pred_mtx));
                
                ref_indices = 1 : (N_ext - N);
                for i = posInVector
                    
                    % optimal weights
                    opt_weights = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, i);
                    
                    % normalize to 1 and compensate other weights
                    opt_weights = (1 - sum(opt_weights))/size(opt_weights, 1) * ones(size(opt_weights, 1), 1) + opt_weights;
                    
                    opt_pred_mtx(i, ref_indices) = opt_weights;
                    
                    weightsRecord(:, count) = opt_weights;
                end
            end
            
            % Visualize weights vs samples
            plot(1 : (N_ext - N), opt_pred_mtx(i, 1 : (N_ext - N)), [lineStyleSet{count}], 'LineWidth', 1.5); hold on;
            legendInfo{count} = ['\sigma = ' num2str(sigma) ];
            h_title = title( ['x_{' num2str(row) ',' num2str(column) '}'] );
            axis([1 17 -0.1 1]);
            h_xlabel = xlabel('reference pixels');
            h_ylabel = ylabel('weights');
            grid on;
            
            set(h_title, 'FontSize', 15, 'FontName','Times New Roman');
            set(h_xlabel, 'FontSize', 15, 'FontName','Times New Roman');
            set(h_ylabel, 'FontSize', 15, 'FontName','Times New Roman');
        end
        
        h_legend = legend(legendInfo);
        set(h_legend, 'FontName','Times New Roman', 'FontSize', 10, ...
            'FontName','Times New Roman');
        set(gca,'FontSize', 12, 'FontName','Times New Roman');
        
    end
end