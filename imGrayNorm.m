function [res] = imGrayNorm(imin, pos_only)
%IMGRAYNORM 
% Normalize input image into 'imdouble' ranges within 0~1.
%
if nargin<2, pos_only = false; end;
%% convert complex pixels into normalized length
comp_mod = imag(imin);
comp_idx = find(comp_mod>0);
if(~isempty(comp_idx))
    imin(comp_idx) = norm(imin(comp_idx));
end
imin = real(imin);  imin = im2double(imin);
%% scalar normalization
gmin = min(imin,[],'all');
gmax = max(imin,[],'all');
if ~pos_only
    gdel = gmax - gmin + eps;
    res = (imin - gmin)./gdel;
    %fprintf('input image: min/max = %g/%g\n', gmin, gmax);
else
    res = exp(imin)./exp(gmax);
end
end

