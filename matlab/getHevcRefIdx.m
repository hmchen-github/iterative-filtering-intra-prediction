% This function is to get the index of pixels for HEVC intra mode.
%
% Input:
%     width: size of the block, e.g., 4 for 4x4 block, 8 for 8x8 block.
%     predIntraMode: intra prediction mode in HEVC (34 mode).
% Output:
%     ref_indices: 1-D vector of index of pixels for |predIntraMode|.
%
function ref_indices = getHevcRefIdx( width, predIntraMode )

if ( predIntraMode == 1)
    ref_indices = [ 1: width + 1, 2 * width + 2 : 3 * width + 1 ];
    return;
end

if ( predIntraMode <= 10)
    ref_indices = 1 : 3 * width + 1;
elseif ( predIntraMode >= 26)
    ref_indices = [ 1 : width + 1, 2 * width + 2 : 4 * width + 1 ];
else
    ref_indices = 1 : 4 * width + 1;
end
end