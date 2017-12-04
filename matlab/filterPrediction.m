function [filter_pred_mtx] = filterPrediction( block_size, pred_mtx, iteration, kernel )

N_ext = block_size^2 + 4 * block_size + 1;

filter_mtx = zeros( N_ext, N_ext );

% ====== generate mapping from block to vector ======
map_size = 2 * block_size + 1;
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
for x = 2 : block_size + 1
    for y = 2 : block_size + 1
        count = count + 1;
        order_map( y, x ) = count;
    end
end


% ====== assign prediction weights ======

for i = 1 : 1 + 4 * block_size
    filter_mtx(i, i) = 1;
end

kernel_size = size(kernel, 1);

if ( sum(sum((abs(kernel)))) == 0 )
    kernel = ones(kernel_size, kernel_size);
end

for y = 2 : block_size + 1
    for x = 2 : block_size + 1
        sum_kernel = 0;
        for off_y = - floor(kernel_size/2) : floor(kernel_size/2)
            for off_x = - floor(kernel_size/2) : floor(kernel_size/2)
                if (y + off_y < 1 || y + off_y > block_size + 1 || x + off_x < 1 || x + off_x > block_size + 1)
                    continue;
                end
                sum_kernel = sum_kernel + kernel(off_y + floor(kernel_size/2) + 1, off_x + floor(kernel_size/2) + 1);
            end
        end
        
        for off_y = - floor(kernel_size/2) : floor(kernel_size/2)
            for off_x = - floor(kernel_size/2) : floor(kernel_size/2)
                if (y + off_y < 1 || y + off_y > block_size + 1 || x + off_x < 1 || x + off_x > block_size + 1)
                    continue;
                end
                filter_mtx( order_map(y, x), order_map(y + off_y, x + off_x) ) = ...
                    kernel(off_y + floor(kernel_size/2) + 1, off_x + floor(kernel_size/2) + 1)/sum_kernel;
            end
        end
    end
end

% ====== apply the filter ======

filter_pred_mtx = filter_mtx^iteration * pred_mtx;