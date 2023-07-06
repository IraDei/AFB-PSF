function [x] = AFB_PSF(img, varargin)
%AFB_PSF 
% Expansion Anisotropic Filter Bank as pixel level features for IRST components. 
% Please refer to this paper for source:
% E. Zhao, W. Zheng, M. Li, Z. Niu, X. Liu and J. Wang, "A Fast Detection
% Method Using Anisotropic Guidance for Infrared Small Target Under Complex
% Scenes," in IEEE Geoscience and Remote Sensing Letters, vol. 20, pp. 1-5,
% 2023, Art no. 6002805, doi: 10.1109/LGRS.2023.3243469.
%  
dbg = 01;
if(ischar(img))
    img = imread(img);
end
[imgR ,imgC, dimension]=size(img);
if(dimension>2) 
    img = rgb2gray(img);
    dimension = 1;
end
img = im2double(img);

vmap = parse_varargin('afb_', varargin{:});
if vmap.isKey('measR')
    measR = vmap('measR');  radius = measR{1};
else
    % default mask radius rank set to 1
    radius = [1;]; 
end

%% Generate PSF filter masks
% Comment: 
% The rotation here seems to be nonsense since the value difference is too
% minute
[psf, psfr] = applyRotatedGaussianPSF(1, 2, 45);
corr_psf = [1 0 1;0 0 0;1 0 1;] .* psfr .* (1./psf) + [0 1 0; 1 1 1; 0 1 0;];

%% Generate Anisotropic Filter Banks
% The proposed AFB constructs four direction filtering kernels, each
% satisfying the characteristics of the edge filtering operator, that is,
% the sum of filtering coefficients is 0.
asf_Bases(:,:,1) = [0 0 0; 1 2 1; 0 0 0;];  % base for 0 degree kernel
asf_Bases(:,:,2) = [0 -1 0; -1 2 -1; 0 -1 0;];  % base for 45 degree kernel
asf_Bases(:,:,3) = [0 1 0; 0 2 0; 0 1 0;];  % base for 90 degree kernel
asf_Bases(:,:,4) = [0 -1 0; -1 2 -1; 0 -1 0;];  % base for 135 degree kernel
% deviation positions
asf_kpos{1} = [1,1; 1,3; 3,1; 3,3;];
asf_kpos{2} = [1,3; 3,1;];
asf_kpos{3} = [1,1; 1,3; 3,1; 3,3;];
asf_kpos{4} = [1,1; 3,3;];

[asf_maskR, asf_maskC, asf_maskN] = size(asf_Bases);
asf_mask = zeros(asf_maskR, asf_maskC, asf_maskN);
afb = zeros(asf_maskR, asf_maskC, asf_maskN);   % Anisotropic filter based PSF
% Directional features by simple conv filtering
df = zeros(imgR, imgC, asf_maskN);
% padding input image for AF conv
imgPad = padarray(img, [1, 1], 0, 'both');
for i = 1:asf_maskN
    asf_mask(:,:,i) = GenAFMask(asf_Bases(:,:,i), asf_kpos{i});
    afb(:,:,i) = asf_mask(:,:,i) .* corr_psf;
    % Conv the afb alone the img, during which the mask afb is always
    % expands its centroid at the indicator position of the object signal.
    %df(:,:,i) = conv2(img, afb(:,:,i), 'same');
    df(:,:,i) = conv2(imgPad, afb(:,:,i), 'valid');
end

% Channel minimize pooling via 1*1 conv
df_min = zeros(imgR, imgC);
de = zeros(imgR, imgC); % directional energy
for y = 1:imgR
    for x = 1:imgC
        df_min(y, x) = min(df(y,x,:));
        de(y, x) = max(sum(df(y,x,:)), 0);
    end
end
% normalization
df_min = imGrayNorm(df_min);
de = imGrayNorm(de);
%multiImgShow({df_min, de});
%% Compute the multi-scale Isotropic measures
radiusN = numel(radius);
im = zeros(imgR, imgC, radiusN);
m = zeros(imgR, imgC, radiusN);
imde = zeros(imgR, imgC, radiusN);
for i = 1:radiusN
    [im(:,:,i), m(:,:,i)] = rangedIM(df_min, radius(i));
    % fuse scaled IMDE saliency
    imde(:,:,i) = imGrayNorm(im(:,:,i) .* de);
end

%% Fuse the final IMDE via 1*1 conv
IMDE_fuse = zeros(imgR, imgC);
for y = 1:imgR
    for x = 1:imgC
        IMDE_fuse(y, x) = max(imde(y,x,:));
    end
end

%% feature package
if dbg
    fhd = multiImgShow({img, df_min, de, im, m, IMDE_fuse});
    set(fhd, 'Name', sprintf('AFB-PSF: radius = %s', mat2str(radius)));
end
x = cat(3, df, df_min, de, im, m, IMDE_fuse);

end

function [kernel, rotatedKernel] = applyRotatedGaussianPSF(frank, sigma, angle, varargin)
    debugDisp = 0;
    % Generate a 2D Gaussian filter kernel
    kernelSize = frank * sigma;  % Adjust the size according to the desired filter strength
    kernelRadius = floor(kernelSize / 2);
    % The filter size equals to frank*2 + 1
    [x, y] = meshgrid(-kernelRadius:kernelRadius, -kernelRadius:kernelRadius);
    kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    kernel = kernel / sum(kernel(:));  % Normalize the kernel
    
    % Rotate the kernel by the specified angle
    rotatedKernel = imrotate(kernel, angle, 'bicubic', 'crop');
    
    % Apply the filter to the image using convolution
    %filteredImage = conv2(image, rotatedKernel, 'same');
    
    % Display filter shape
    if debugDisp
        fgd = figure; 
        localSubPlot(fgd, 1, kernel, 'Gaussian Filter: PSF');
        localSubPlot(fgd, 2, rotatedKernel, 'Gaussian Filter: PSF-45');
    end
    function localSubPlot(fig, idx, img, ftitle)
        figure(fig);
        subplot(1,2,idx), surf(x, y, img); colormap('jet');
        shading interp; colorbar;
        title(ftitle);
        xlabel('X');    ylabel('Y');        zlabel('Magnitude');
    end
end

function asf = GenAFMask(base, kpos, varargin)
% Generate Anisotropic filter mask by filling compensation items into given
% preset base mask according to positions given in kpos.
% The grid positions are suggested in matrix coords.
% The base mask shall not contain any negative values.
%

%[mR, mC] = size(base);  
kposN = size(kpos, 1);  % compensation grid positions in rows of [row, col]
asf = base;
if kposN >0 
    kvSum = sum(base, 'all');   kval = - kvSum / kposN;
    for i = 1:kposN
        asfR = kpos(i, 1);  asfC = kpos(i, 2);
        asf(asfR, asfC) = kval;
    end
end
end

function [im, m] = rangedIM(df, r, varargin)
% Compute the isotropic measure under given radius.
% The filter size equals to 2*r+1.
%
maskL = r*2 + 1;
accMask = ones(maskL, maskL);
% compute local block sum for each position
dfPad = padarray(df, [r, r], 0, 'both');
dfRSum = conv2(dfPad, accMask, 'valid');

m = df ./ dfRSum;
% Compute the single point isotropic measure
imPnt = -m .* log(1 - m.^2);
imPntPad = padarray(imPnt, [r, r], 0, 'both');
im = conv2(imPntPad, accMask, 'valid');

%multiImgShow({dfRSum, m, imPnt, im});
end
