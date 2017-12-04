% This function is to get the optimal weights for the prediction in a
% mse sense. This method is based on Jayant's "Digital Coding of Waveforms"
% page 269
%
% Input:
%     cov_mtx_ext: covariance matrix (assuming zero-mean)
%     ref_indices: indices of the references points
%     pred_index: the index of the to-be-predicted point
% Output:
%     opt_weights: optimal weights
%


function [opt_weights] = getOptimalPredictionWeights(cov_mtx_ext, ref_indices, pred_index)

rxx = cov_mtx_ext(ref_indices, pred_index);

Rxx = cov_mtx_ext(ref_indices, ref_indices);

opt_weights = inv(Rxx) * rxx;

end