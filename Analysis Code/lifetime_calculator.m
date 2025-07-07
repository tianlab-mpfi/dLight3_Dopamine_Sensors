%function allows for lifetime calculation of a single frame

function lt = lifetime_calculator(arr, total_time_window, prctile_used)
    dim = size(arr);
    pixels = dim(1);
    time_bins_available = dim(3) - 1;
    time_step = total_time_window/dim(3);
    time_arr = 0:time_step:(dim(3)-3)*time_step;
    subimg = arr(:,:,1:time_bins_available); %take the first 31 elements since the last element of array is sum
    prctile_33 = prctile(arr(:,:,end), prctile_used, 'all'); %find the 33rd percentile (exclude lower pixels)
    [i, j] = find(arr(:,:,end) < prctile_33); %get the pixels (in 2d) that are less than 33rd percentile, top 66% of pixels (Yao's paper)
    elements_to_zero = [i,j]; %make an array out of these pixel indices
    ind_to_zero = sub2ind(size(arr(:,:,end)), elements_to_zero(:,1), elements_to_zero(:,2)); %convert the 2d array of indices to 1d index array
    
    arr_1 = ones(pixels,pixels);
    arr_1(ind_to_zero) = 0; %create a 2d mask of elements you want to include
    arr_2 = repmat(arr_1, [1,1,time_bins_available]); %make this into 3d based on time bins available
    zerod_subimg = subimg .* double(arr_2); %mask out the lowest pixel values
    
    pooled_zerod_subimg = sum(zerod_subimg, [1,2]); %pool all pixels together in image
    ds_pooled_zerod_subimg = pooled_zerod_subimg(1,:);
    lt = sum(ds_pooled_zerod_subimg(2:time_bins_available).*time_arr)/sum(ds_pooled_zerod_subimg(2:time_bins_available)); %exclude first element before the peak and calculate average lifetime
    
    disp(lt);%display average lifetime