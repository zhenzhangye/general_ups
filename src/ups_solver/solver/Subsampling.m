%% Subsample the images if needed
%  input: data.I            images             h*w*c*N
%         data.mask         binary             h*w
%         data.K            camera intrinsics  3*3
%         options.ratio     subsample ratio    int
%  output: all data after subsampling.
function [data] = Subsampling(data, options)

if(options.ratio>1)
    data.I        = data.I(1:options.ratio:end,1:options.ratio:end,:,:);
    
    data.mask     = data.mask(1:options.ratio:end,1:options.ratio:end);
 
    data.K(1:2,:) = data.K(1:2,:)./options.ratio;
end

end