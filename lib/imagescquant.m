% Plots a grayscale image with intensity axis scaled to the quantile for
% cumulative probability specified in quant.
function imagescquant(img,quant)
    if nargin < 2
        quant = 0.95;
    end
    imagesc(img); axis equal; axis off; colormap gray;
    imlims = [0,quantile(as_vector(img),quant)];
    caxis(imlims);
end