% This function is to get the prediction matrix and covariance matrix of
% the signal of the original block.

function [pred_mtx, cov_mtx_ext] = getHevcIntraPredAndExtCovMtx( BlockSize, ...
                                                                 predModeIntra, ...
                                                                 rho, ...
                                                                 alpha, ...
                                                                 eta, ...
                                                                 sigma, ...
                                                                 disableIntraBoundaryFilter...
                                                                )

intraPredModeSet = {'INTRA_PLANAR',...
    'INTRA_DC',...
    'INTRA_ANGULAR2',...
    'INTRA_ANGULAR3',...
    'INTRA_ANGULAR4',...
    'INTRA_ANGULAR5',...
    'INTRA_ANGULAR6',...
    'INTRA_ANGULAR7',...
    'INTRA_ANGULAR8',...
    'INTRA_ANGULAR9',...
    'INTRA_ANGULAR10',...
    'INTRA_ANGULAR11',...
    'INTRA_ANGULAR12',...
    'INTRA_ANGULAR13',...
    'INTRA_ANGULAR14',...
    'INTRA_ANGULAR15',...
    'INTRA_ANGULAR16',...
    'INTRA_ANGULAR17',...
    'INTRA_ANGULAR18',...
    'INTRA_ANGULAR19',...
    'INTRA_ANGULAR20',...
    'INTRA_ANGULAR21',...
    'INTRA_ANGULAR22',...
    'INTRA_ANGULAR23',...
    'INTRA_ANGULAR24',...
    'INTRA_ANGULAR25',...
    'INTRA_ANGULAR26',...
    'INTRA_ANGULAR27',...
    'INTRA_ANGULAR28',...
    'INTRA_ANGULAR29',...
    'INTRA_ANGULAR30',...
    'INTRA_ANGULAR31',...
    'INTRA_ANGULAR32',...
    'INTRA_ANGULAR33',...
    'INTRA_ANGULAR34'};

intraPredAngleSet = [NaN, NaN, ...% Planar, DC
     32,  26,  21,  17,  13,  9,  5,  2, 0, -2, -5, -9, -13, -17, -21, -26, ...   %INTRA_ANGULAR2 ~ INTRA_ANGULAR17
    -32, -26, -21, -17, -13, -9, -5, -2, 0,  2,  5,  9,  13,  17,  21,  26, 32]; %INTRA_ANGULAR18 ~ INTRA_ANGULAR34

invAngleSet = round( 2^13./intraPredAngleSet );
    
nS = BlockSize;

%% Check
if (predModeIntra < 0 || predModeIntra > 34)
    disp 'ATTENTION : mode_index must be between 0 to 34';
    pred_mtx = 0;
    cov_mtx_ext = 0;
    return;
end

%% Step 1
% Generate a big vector including block pixels, p(x,-1) and p(-1,x), x=0..nS-1.
% record corresponding pixel's loaction [x,y]

p_vec = [];
vec_num = 0;
%********************
%for i = -1:(-1+2*nS)
%    VecNum = VecNum + 1;
%    pVec = [pVec; [-1 i]];
%end
%********************

% Left
for y = -1 : ( -1 + 2 * nS)
    vec_num = vec_num + 1;
    p_vec = [p_vec; [-1 y]];
end

% Top
for x = 0 : ( -1 + 2 * nS)
    vec_num = vec_num + 1;
    p_vec = [p_vec; [x -1]];
end

% Main body
for x = 0 : ( -1 + nS)
    for y = 0 : ( -1 + nS)
        vec_num = vec_num + 1;
        p_vec = [p_vec; [x y]];
    end
end

%vec_num
%p_vec

%% Step 2
% Generate the covariance matrix of the big vector
cov_mtx_ext = zeros( vec_num, vec_num );

for i = 1 : vec_num
    for j = 1 : vec_num
        ix = p_vec(i, 1); iy = p_vec(i, 2);
        jx = p_vec(j, 1); jy = p_vec(j, 2);
        cov_mtx_ext(i,j) = ...
            rho^sqrt(((ix - jx)*cos(alpha) - (iy - jy)*sin(alpha))^2 + ...
            eta^2*((ix - jx)*sin(alpha) + (iy - jy)*cos(alpha))^2);
    end
end

for i = 1 : 4 * nS + 1
    cov_mtx_ext(i, i) = cov_mtx_ext(i, i) + sigma^2;
end


%% ====== generate mapping from block to vector ======
map_size = 2 * nS + 1;
order_map = zeros( map_size, map_size );
count  = 0;
for y = 1 : map_size
    count = count + 1;
    order_map(y, 1) = count;
end
for x = 2 : map_size
    count = count + 1;
    order_map(1, x) = count;
end
for x = 2 : nS + 1
    for y = 2 : nS + 1
        count = count + 1;
        order_map( y, x ) = count;
    end
end

%% Step 3
% Generate the intraPrediction maxtrix pred_mtx that
% pred_mtx*pVec = pVec - pred_pVec

pred_mtx = zeros( vec_num, vec_num );

for i = 1 : 1 + 4 * nS
    pred_mtx(i, i) = 1;
end

% VER Mode
if (predModeIntra == 26)
    for i = (4*nS + 2) : vec_num
        ix = p_vec(i, 1);
        pred_pos_in_p_vec = order_map(1, ix + 2);
        pred_mtx(i, pred_pos_in_p_vec) = 1;
    end
    
    if ( disableIntraBoundaryFilter || nS > 32 )
        return;
    else
        % filtering along boundaries
        % [0][y]
        x = 2;
        for y = 2 : nS + 1
            pred_mtx( order_map(y, x), :) = 0; % clear previous weights
            pred_mtx( order_map(y, x), order_map(1, x) ) = 1;
            pred_mtx( order_map(y, x), order_map(y, 1) ) = 1/2;
            pred_mtx( order_map(y, x), order_map(1, 1) ) = -1/2;
        end
        return;
    end
end


% HOR Mode
if (predModeIntra == 10)
    for i = (4*nS + 2) : vec_num
        iy = p_vec(i,2);
        pred_pos_in_p_vec = order_map(iy + 2, 1);
        pred_mtx(i, pred_pos_in_p_vec) = 1;
    end
    
    if ( disableIntraBoundaryFilter || nS > 32 )
        return;
    else
        % filtering along boundaries
        % [x][0]
        y = 2;
        for x = 2 : nS + 1
            pred_mtx( order_map(y, x), :) = 0; % clear previous weights
            pred_mtx( order_map(y, x), order_map(y, 1) ) = 1;
            pred_mtx( order_map(y, x), order_map(1, x) ) = 1/2;
            pred_mtx( order_map(y, x), order_map(1, 1) ) = -1/2;
        end

        return;
    end
end

% DC Mode
% Note that a filter after DC prediction is not applied!
if (predModeIntra == 1)
    k = log2(nS);
    for i = (4*nS + 2) : vec_num
        for pred_pos_in_p_vec = 2 : nS+1 % top
            pred_mtx(i, pred_pos_in_p_vec) = 1/(2^(k+1));
        end
        for pred_pos_in_p_vec = 2*nS+2 : 3*nS+1 % left
            pred_mtx(i, pred_pos_in_p_vec) = 1/(2^(k+1));
        end
    end
    
    if ( disableIntraBoundaryFilter || nS > 32 )
        return;
    else
        % filtering along boundaries
        boundary_filter_mtx = zeros(vec_num, vec_num);
        % [-1][y] & [x][-1]
        boundary_filter_mtx(1 : 4 * nS + 1, 1 : 4 * nS + 1) = eye(4 * nS + 1);
        % [0][0]
        x = 2;
        y = 2;
        boundary_filter_mtx( order_map(y, x), order_map(y - 1, x) ) = 1/4;
        boundary_filter_mtx( order_map(y, x), order_map(y, x) ) = 2/4;
        boundary_filter_mtx( order_map(y, x), order_map(y, x - 1) ) = 1/4;
        % [x][0]
        y = 2;
        for x = 3 : nS + 1
            boundary_filter_mtx( order_map(y, x), order_map(y, x) ) = 3/4;
            boundary_filter_mtx( order_map(y, x), order_map(y - 1, x) ) = 1/4;
        end
        % [0][y]
        x = 2;
        for y = 3 : nS + 1
            boundary_filter_mtx( order_map(y, x), order_map(y, x) ) = 3/4;
            boundary_filter_mtx( order_map(y, x), order_map(y, x - 1) ) = 1/4;
        end
        
        for x = 3 : nS + 1
            for y = 3 : nS + 1
                boundary_filter_mtx( order_map(y, x), order_map(y, x) ) = 1;
            end
        end
        
        pred_mtx = boundary_filter_mtx * pred_mtx;
        return;
    end
end


%% Step 4 (Angular Mode)
% if mode_index >= 3

%Step 4.1
% Derive refMain[x] (x = -nS...2*nS)specifing the reference samples' location
% We need a Offset=nS+1 here so that the x+Offset are always > 0


Offset = nS+1;
refMain = zeros(2*nS+Offset, 2); %refMain(:,1)--(x,:), refMain(:,2)--(:,y)

intraPredAngle = intraPredAngleSet( predModeIntra + 1 );
invAngle = invAngleSet( predModeIntra + 1 );

if predModeIntra >= 18
    for refIndex = 0 : nS
        refMain(refIndex+Offset, :) = [-1 + refIndex, -1];
    end
    if intraPredAngle < 0
        for refIndex = floor((nS*intraPredAngle)/32) : -1
            refMain(refIndex+Offset, :) = [-1, -1+floor((refIndex*invAngle+128)/256)];
        end
    else
        for refIndex = nS+1 : 2*nS
            refMain(refIndex+Offset, :) = [-1 + refIndex, -1];
        end
    end
else   
    for refIndex = 0 : nS
        refMain(refIndex+Offset,:) = [-1, -1 + refIndex];
    end
    if intraPredAngle < 0
        for refIndex = floor((nS*intraPredAngle)/32) : -1
            refMain(refIndex+Offset, :) = [-1+floor((refIndex*invAngle+128)/256), -1];
        end
    else
        for refIndex = nS+1:2*nS
            refMain(refIndex+Offset, :) = [-1, -1 + refIndex];
        end
    end
end

% Step 4.2
%****************************************************************
% iIdx = ((y+1)*intraPredAngle) >> 5
% iFact = ((y+1)*intraPredAngle) & 31
% if iFact <> 0
%   predSample[x,y] = ((32-iFact)*refMain[x+iIdx+1] + iFact*refMain[x+iIdx+2] + 16) >> 5
% else
%   predSample[x,y] = refMain[x+iIdx+1]
% end
%****************************************************************

for i = (4*nS + 2) : vec_num
    if predModeIntra >= 18
        x = p_vec(i,1);
        y = p_vec(i,2);
    else
        x = p_vec(i,2);
        y = p_vec(i,1);
    end
    iIdx = floor((y+1)*intraPredAngle/32);
    if (y+1)*intraPredAngle >= 0
        iFact = bitand((y+1)*intraPredAngle, 31);
    else
        iFact = bitand(bitxor(-(y+1)*intraPredAngle,2^20-1)+1, 31);
    end
    %iFact
    if iFact ~= 0
        if refMain(x + iIdx + 1 + Offset, 1) == -1
            pred_pos_in_p_vec_1 = refMain(x + iIdx + 1 + Offset, 2) + 2;
        else
            pred_pos_in_p_vec_1 = 2*nS + 2 + refMain(x + iIdx + 1 + Offset, 1);
        end
        
        if refMain(x+iIdx+2+Offset, 1) == -1
            pred_pos_in_p_vec_2 = refMain(x + iIdx + 2 + Offset, 2) + 2;
        else
            pred_pos_in_p_vec_2 = 2*nS + 2 + refMain(x + iIdx + 2 + Offset, 1);
        end
        pred_mtx(i, pred_pos_in_p_vec_1) = (32 - iFact)/32;
        pred_mtx(i, pred_pos_in_p_vec_2) = iFact/32;
    else
        if refMain(x + iIdx + 1 + Offset, 1) == -1
            pred_pos_in_p_vec = refMain(x + iIdx + 1 + Offset, 2) + 2;
        else
            pred_pos_in_p_vec = 2*nS + 2 + refMain(x + iIdx + 1 + Offset,1);
        end
        pred_mtx(i, pred_pos_in_p_vec) = 1;
    end
end
