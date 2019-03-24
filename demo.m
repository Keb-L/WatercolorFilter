function varargout = demo(varargin)
% DEMO M-file for demo.fig
%      DEMO, by itself, creates a new DEMO or raises the existing
%      singleton*.
%
%      H = DEMO returns the handle to a new DEMO or the handle to
%      the existing singleton*.
%
%      DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO.M with the given input arguments.
%
%      DEMO('Property','Value',...) creates a new DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help demo

% Last Modified by GUIDE v2.5 03-Dec-2017 11:03:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% declare global variables to hold the image and handle
global X;
global hAxes1;
global hAxes2;

% open an image
[FileName,PathName] = uigetfile('*.bmp;*.tif;*.jpg;*.hdf','Select the image file');
FullPathName = [PathName,'\',FileName];
X = imread(FullPathName);

% get the handle and display it
hAxes1 = findobj(gcf,'Tag','axes1');
hAxes2 = findobj(gcf,'Tag','axes2');

%set(gcf, 'CurrentAxes', hAxes1);
axes(handles.axes1);
imshow(X);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;
global hAxes1;
global hAxes2;

hAxes1 = findobj(gcf,'Tag','axes2');
hAxes2 = findobj(gcf,'Tag','axes1');

%compile_edison_wrapper
%figure;
%imshow(I);

nrow = size(X, 1);      % number of rows
ncol = size(X, 2);      % number of columns

thresh = 0.0004;        % percentage of pixel space that a segmentation must use
beta = 1.5;               % global darkening constant for texture
beta_edge = 0.5;        % darkening constant for edghes
%sigma = 2;              % gaussian kernel sigma
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

%figure;
%imshow(texture);
tic
disp('Initialized constants...');

% EDISON parameters
% segments image
minArea = floor(thresh * nrow * ncol);
[fimg labels modes regsize grad conf] = edison_wrapper(X,@rgb2lab, 'MinimumRegionArea', minArea,...
                                            'RangeBandWidth', 2,... 
                                            'SpatialBandWidth', 11,... 
                                            'SpeedUp', 3,...
                                            'GradientWindowRadius', 3,...
                                            'MixtureParameter', 0.3,...
                                            'EdgeStrengthThreshold', 0.3);


disp('Segmentation complete...');
img_seg = lab2rgb(fimg);
%figure;
%hAxes1 = findobj(gcf,'Tag','axes1');
%set(gcf, 'CurrentAxes', hAxes1);
%imshow(img_seg);

num_region = length(regsize);
%segmented_images = cell(1, num_region);
%segmented_binary = cell(1, num_region);
%segmented_morph = cell(1, num_region);

% utilizes morphological smoothing to smooth edges of each segment
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

disp('Smoothing complete...');
%figure;
%imshow(img_morph);

% Finds edges
H = edgeKernel(edge_N);
img_V = conv3d(img_morph, H);
img_H = conv3d(img_morph, H');
img_mag = rgb2gray(sqrt(img_H.^2 + img_V.^2));

img_mag(img_mag<edge_thresh) = 0;

%figure;
%imshow(img_mag);

% Adds turbulence
img = edgeDarken(img_morph, texture, beta);
disp('Turbulence added...');
% Darkens strong edges
img = edgeDarken(img, img_mag, beta_edge);
disp('Edges modified...');

% Reduces saturation of the image
img = gamma(img, gamma1, gamma2);
disp('Gamma correction complete...');
toc

%set(gcf, 'CurrentAxes', hAxes2);
axes(handles.axes2);
imshow(img);


function im = turbulence(sz, iter)
    noise = randn(sz);
    im = noise / iter;
   
    while iter >= 1
        noise = noise(1:ceil(end/2), 1:ceil(end/2));
        im = im + imresize(noise, sz) ./ iter; 
        
        iter = iter/2;    
    end
    %normalize values
     im = (im - min(min(im(:,:)))) ./ (max(max(im(:,:))) - min(min(im(:,:))));

% applies gamma correction in the saturation and value domain
% power law transform saturation by gamma1, value by gamma2
function im = gamma(img, gamma1, gamma2)
img = rgb2hsv(img);
img_mask1 = img(:,:,1);
img_mask2 = img(:,:,2).^gamma1;
img_mask3 = img(:,:,3).^gamma2;
im(:,:,1) = img_mask1;
im(:,:,2) = img_mask2;
im(:,:,3) = img_mask3;
im = hsv2rgb(im);

% Performs a convolution on a 3D matrix with a specified 2D Kernel
function im = conv3d(img, H)
im = img;
img_mask1 = conv2(img(:,:,1), H, 'same');
img_mask2 = conv2(img(:,:,2), H, 'same');
img_mask3 = conv2(img(:,:,3), H, 'same');
im(:,:,1) = img_mask1;
im(:,:,2) = img_mask2;
im(:,:,3) = img_mask3;

% Applies pigment modification as described in Bousseau's paper
function im = edgeDarken(img, I, beta)
d = 1 + beta * (I - 0.5);
im = im2double(img);
im = im - (im - im.^2).*(d-1);

% Generates the edge gradient kernel as described in Doran's paper
function H = edgeKernel(N)
H = zeros(N, N);
for i = 1:floor(N/2)
    H(i, :) = [i:1:ceil(N/2)+i-1 ceil(N/2)-2+i:-1:i];
end
H = -H + flipud(H);
