function [fhd] = multiImgShow(imCell, varargin)
%MULTIIMGSHOW 
% Display images in 'imCell' in a single plot with automatical subplot
% sizes.
%   

fhd = figure;
imQty = numel(imCell);  subpC = ceil(sqrt(imQty));  subpR = floor(imQty/subpC) + 1;
for i = 1:imQty
    subplot(subpR, subpC, i), imshow(imCell{i},[]);
end

end

