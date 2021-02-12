
%% Skeleton code: Denoising IHS for CT reconstruction with Gaussian noise
%|
%| A. Perelli, M. S. Andersen, 
%| "Regularization by Denoising Sub-sampled Newton Method for 
%| Spectral CT Multi-Material Decomposition" 
%| Philosophical Transactions of the Royal Society A, (2021)
%| DOI: 

%% ********************  Description  *************************
%| This script reconstruct 2D images form mono energetic CT measurements
%| with Gaussian noise.

%% ****************  Software Requirement  ********************
%| For reconstruction, the code requires the ASTRA Toolbox:
%| https://github.com/astra-toolbox/astra-toolbox
%|
%| W. van Aarle, W. J. Palenstijn, J. Cant, E. Janssens, 
%| F. Bleichrodt, A. Dabravolski, J. De Beenhouwer, 
%| K. J. Batenburg, and J. Sijbers, 
%| "Fast and Flexible X-ray Tomography Using the ASTRA Toolbox" 
%| Optics Express, 24(22), 25129-25147, (2016) 
%| http://dx.doi.org/10.1364/OE.24.025129

%% ************************************************************
%| Copyright 08-02-2021, Alessandro Perelli 
%| Technical University of Denmark (DTU)

clear, close; clc;

% Astra Toolbox
if ispc
    addpath(genpath('astra_toolbox_matlab_win'));
elseif isunix
    addpath(genpath('../../../astra'))
end


% File must be in active directory
path = './test_images/';
filename = 'starfish.tif';

% Spot Toolbox
addpath spot;
load lev_scores_A.mat

addpath(genpath('Denoisers'));

%% Parameters
% Image dimension
d  = 64;
dv = d^2;       % Vectorized dimensions

%% Generate High and Low resolution image forward CT operators
vol_geom = astra_create_vol_geom(d, d);

anglesNumb    = 180;                          % number of projection angles
anglesDegrees = linspace(0, 180, anglesNumb); % projection angles
Detectors     = round(sqrt(2)*d);            % number of detectors
n             = anglesNumb * Detectors;

% using ASTRA-toolbox to set the projection geometry (parallel beam)
proj_geom = astra_create_proj_geom('parallel', 1, Detectors, anglesDegrees*pi/180);

% projection 'line' model used
proj_id   = astra_create_projector('line', proj_geom, vol_geom);

% %% Generate Measurements (sinogram)
% % Input high resolution image - Ground truth
% I = 255*phantom(d);

%% read the original image
fprintf('Reading %s image...', filename);
orig_im = imread(['./test_images/'  filename]);
I = double(orig_im(1:4:end,1:4:end,1));

%% Generate projection data
% Create operator for ASTRA using the GPU.
A = @(x) astra_create_sino_handle(x, proj_id);

% Linear measurement z = Ax
sinogram = A(I);

% % adding Poisson noise
% dose =  1e4;                                         % photon flux (controls noise level)
% sinoNoise = astra_add_noise_to_sino(sinogram, dose); % adding Poisson noise

% add Gaussian noise
variance  = .1;
sinoNoise = sinogram + sqrt(variance)*randn(size(sinogram));

%% Sketching Astra Operators
% Sketched Forward Operator
SA  = @(x, proj_id_SUB) astra_create_sino_handle(x, proj_id_SUB);
% Sketched Backward Operator
SAT = @(y, proj_id_SUB) astra_create_backprojection_handle(y, proj_id_SUB);

% Full Backward Operator
AT = @(y) astra_create_backprojection_handle(y, proj_id);

%% Leverage Scores
m = round(dv/Detectors);
block_lev = zeros(anglesNumb, 1);
for i = 1:anglesNumb
    block_lev(i) = sum( leveragescores( (i-1)*Detectors+1:i*Detectors ) );
end
prob_lev = block_lev / sum(block_lev);

% MSE error function
f_mse  = @(x) mean( I(:) - x(:) ).^2; 

% x^{t+1}  = arg min_x [1/2(x - x^t)^T H^t(x^t)(x - x^t) + < dg(x^t), x - x^t >]
% H^t(x^t) = (S^t[d^2 f(x^t)]^{1/2})^T S^t[d^2 f(x^t)]^{1/2} + d^2 r(x^t)
% \/ g(x) = A'(Ax - b) + 1/v * (D(x) - x)
lambda = .1;

% Denoising function
denoiser =  'DWT'; % 'BM3D'; % 'TV'; % U-Net
denoi = @(noisy) denoise(noisy,sqrt(variance), d, d, denoiser);    

grad_f = @(x, x_denoi) AT( A(x) - sinoNoise ) + lambda*( x_denoi - x );
% hess_f = @(x_denoi, x, proj_id_SUB) SAT(SA(x, proj_id_SUB)) + lambda*( div_denoi(x_denoi,x,lambda,denoi) - 1 );

%% Iterative Hessian Sketching
% Initialization
max_iter_IHS = 100;
max_iters_cg = 100;
mse_ihs      = zeros(max_iter_IHS, 1);
x            = zeros(d);

for i = 1:max_iter_IHS
    
    % Random Leverage Scores Sketching
    sub_idx       = randsample(anglesNumb, m, true, prob_lev);
    sub_idx       = sort(unique(sub_idx));           % eliminate duplicates

    proj_geom_SUB = proj_geom;
    proj_geom_SUB.ProjectionAngles = anglesDegrees(sub_idx);

    proj_id_SUB = astra_create_projector('line', proj_geom_SUB, vol_geom);

    % Conjugate Gradient
    mse_cgs = zeros(max_iters_cg, 1);

    % Denoise
    x_denoi = denoi(x);
    
    g = zeros(d, 1);
    r = -grad_f(x, x_denoi);

    p = r;

    for k = 1:max_iters_cg

        q     = SA(p, proj_id_SUB);
        delta = norm(q)^2 + lambda*norm(p)^2;
        alpha = norm(r)^2 / delta;
        g     = g + alpha*p;

        gamma = norm(r)^2;

        r     = r + alpha*( SAT(q, proj_id_SUB) + lambda*( div_denoi(x_denoi,x,lambda,denoi) - 1 )*p );
        beta  = norm(r)^2 / gamma;

        p     = r + beta*p;

        % MSE error CG
        mse_cgs(k) = f_mse(g);  

    end
    
    x = x + g;
    
    % MSE error
    mse_ihs(i) = f_mse(x);     
end
