clear all;
close all;

[FileName,PathName] = uigetfile('*.bmp;*.tif;*.jpg;*.hdf','Select the image file');
FullPathName = [PathName,'\',FileName];
I = imread(FullPathName);

%compile_edison_wrapper
figure;
imshow(I);

nrow = size(I, 1);      % number of rows
ncol = size(I, 2);      % number of columns

thresh = 0.0004;        % percentage of pixel space that a segmentation must use [default 0.0004]
beta = 1.5;             % global darkening constant for texture
beta_edge = 0.5;        % darkening constant for edges
sigma = 2;              % gaussian kernel sigma
se_rad = 5;             % structuring element radius
edge_N = 3;             % edge gradient size
gamma1 = 1.5;           % S gamma
gamma2 = 0.8;         	% V gamma
turb_iter = 32;         % turbulence degree, higher value = lower turbulence, more smooth

edge_thresh = 0.3;

%n = ceil(rand(1) * 6);
%texture = im2double(rgb2gray(imread([pwd,'\texture\texture',num2str(n),'.jpg'])));
%texture = imresize(texture, [nrow ncol]);
texture = turbulence([nrow, ncol], turb_iter);
texture = texture.^0.4;

figure;
imshow(texture);

% EDISON parameters
minArea = floor(thresh * nrow * ncol);
tic
[fimg labels modes regsize grad conf] = edison_wrapper(I,@rgb2lab, 'MinimumRegionArea', minArea,...
                                            'RangeBandWidth', 2,... 
                                            'SpatialBandWidth', 11,... 
                                            'SpeedUp', 3,...
                                            'GradientWindowRadius', 3,...
                                            'MixtureParameter', 0.3,...
                                            'EdgeStrengthThreshold', 0.3);

img_seg = lab2rgb(fimg);
figure;
imshow(img_seg);

num_region = length(regsize);
%segmented_images = cell(1, num_region);
%segmented_binary = cell(1, num_region);
%segmented_morph = cell(1, num_region);

SE = strel('disk', se_rad);
img_morph = fimg;
for i = 1:num_region
    %segmented_binary{i} = (labels == i-1);
    segmented_binary = (labels == i-1);
    %segmented_morph{i} = imclose(segmented_binary, SE);
    segmented_morph = imclose(segmented_binary, SE);

    img_mask1 = img_morph(:,:,1);
    img_mask2 = img_morph(:,:,2);
    img_mask3 = img_morph(:,:,3);
    
    img_mask1(segmented_morph) = modes(1, i);
    img_mask2(segmented_morph) = modes(2, i);
    img_mask3(segmented_morph) = modes(3, i);
    
    img_morph(:,:,1) = img_mask1;
    img_morph(:,:,2) = img_mask2;
    img_morph(:,:,3) = img_mask3;

end

img_morph = lab2rgb(img_morph);

figure;
imshow(img_morph);

H = edgeKernel(edge_N);
img_V = conv3d(img_morph, H);
img_H = conv3d(img_morph, H');
img_mag = rgb2gray(sqrt(img_H.^2 + img_V.^2));

img_mag(img_mag<edge_thresh) = 0;

figure;
imshow(img_mag);

img = edgeDarken(img_morph, texture, beta);
img = edgeDarken(img, img_mag, beta_edge);

img = gamma(img, gamma1, gamma2);

figure;
imshow(img);
toc

function im = turbulence(sz, iter)
    noise = randn(sz);
    im = noise / iter;
   
    while iter >= 1
        noise = noise(1:ceil(end/2), 1:ceil(end/2));
        im = im + imresize(noise, sz) ./ iter; 
        
        iter = iter/2;    
    end
    
     im = (im - min(min(im(:,:)))) ./ (max(max(im(:,:))) - min(min(im(:,:))));
end

function im = gamma(img, gamma1, gamma2)
img = rgb2hsv(img);
img_mask1 = img(:,:,1);
img_mask2 = img(:,:,2).^gamma1;
img_mask3 = img(:,:,3).^gamma2;
im(:,:,1) = img_mask1;
im(:,:,2) = img_mask2;
im(:,:,3) = img_mask3;
im = hsv2rgb(im);
end

function im = conv3d(img, H)
im = img;
img_mask1 = conv2(img(:,:,1), H, 'same');
img_mask2 = conv2(img(:,:,2), H, 'same');
img_mask3 = conv2(img(:,:,3), H, 'same');
im(:,:,1) = img_mask1;
im(:,:,2) = img_mask2;
im(:,:,3) = img_mask3;
end

function im = edgeDarken(img, I, beta)
d = 1 + beta * (I - 0.5);
im = im2double(img);
im = im - (im - im.^2).*(d-1);
end

function H = edgeKernel(N)
H = zeros(N, N);
for i = 1:floor(N/2)
    H(i, :) = [i:1:ceil(N/2)+i-1 ceil(N/2)-2+i:-1:i];
end
H = -H + flipud(H);
end