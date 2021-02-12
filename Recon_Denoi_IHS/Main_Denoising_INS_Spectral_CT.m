
%% Denoising Iterative Hessian Sketching Computed Tomography
%|
%| A. Perelli, M. S. Andersen, 
%| "Regularization by Denoising Sub-sampled Newton Method for 
%| Spectral CT Multi-Material Decomposition" 
%| Philosophical Transactions of the Royal Society A, (2021)
%| DOI: 

%% ********************  Description  *************************
%| This code performs the Spectral CT reconstruction using the Densoising
%| Iterative Hessian Sketching algorithm.

%% ****************  Software Requirement  ********************
%| The code requires:
%| - ASTRA Toolbox: https://github.com/astra-toolbox/astra-toolbox
%|   W. van Aarle, W. J. Palenstijn, J. Cant, E. Janssens, 
%|   F. Bleichrodt, A. Dabravolski, J. De Beenhouwer, 
%|   K. J. Batenburg, and J. Sijbers, 
%|   "Fast and Flexible X-ray Tomography Using the ASTRA Toolbox" 
%|   Optics Express, 24(22), 25129-25147, (2016) 
%|   http://dx.doi.org/10.1364/OE.24.025129
%| - TomoPhantom (to generate the material phantom model)
%| - ToMoBAR: https://github.com/dkazanc/multi-channel-X-ray-CT

%% ************************************************************
%| Copyright 08-02-2021, Alessandro Perelli 
%| Technical University of Denmark (DTU)
%  ************************************************************ 

clear, close; clc;

% Astra Toolbox
if ispc
    addpath(genpath('astra_toolbox_matlab_win'));
elseif isunix
    addpath(genpath('../../../astra'))
end

% load generated data: TomoPhantom model n. 4
load(sprintf(['..' filesep 'SpectralDataGeneration' filesep 'SpectralData_noCrime.mat'], 1i)); 

% adding paths
addpath(sprintf(['..' filesep 'SupplementaryPackages' filesep], 1i));
addpath(sprintf(['..' filesep 'SupplementaryPackages' filesep 'PhotonAttenuation' filesep], 1i));
addpath(sprintf(['..' filesep 'SupplementaryPackages' filesep 'spot' filesep], 1i));
addpath(sprintf(['..' filesep 'SupplementaryPackages' filesep 'gendist' filesep], 1i));

% adding paths denoisers
addpath(sprintf(['..' filesep 'SupplementaryPackages' filesep 'BM3D' filesep], 1i));

% Set the geometry
B       = -log(bsxfun(@times, Y+(Y==0), 1./sb'));
[Ut,Vl] = geocore_phantom(n,Em,'phantom');
Ut_vec  = [reshape(Ut(:,:,1),n^2,1), reshape(Ut(:,:,2),n^2,1), reshape(Ut(:,:,3),n^2,1), reshape(Ut(:,:,4),n^2,1)];
Phantom_nbins = Ut_vec*Vl';

theta = (0:p-1)*360/p;   % projection angles
dom_width   = 1.0;       % width of domain in cm
src_to_rotc = 3.0;       % dist. from source to rotation center
src_to_det  = 5.0;       % dist. from source to detector
det_width   = 2.0;       % detector width
vol_geom  = astra_create_vol_geom(n,n);

% Projection geometry
proj_geom = astra_create_proj_geom('fanflat', n*det_width/nd, nd, (pi/180)*theta,...
            n*src_to_rotc/dom_width, n*(src_to_det-src_to_rotc)/dom_width);

% selecting the number of channels to reconstruct (energy window)
start_channel = 1; end_channel = size(Y,2);
sino = B(:,start_channel:end_channel); 
num_channels = size(sino,2);

% weights for the PWLS model
W = Y/max(Y(:)); 
% Data normalisation 
Yt = Y + (Y==0);
snr_ray = max(1e-5,-log(bsxfun(@rdivide,Yt,sb'))).*sqrt(Yt);


%% Generate Forward operator
% Create operator for ASTRA using the GPU.
A = @(x) astra_create_sino_handle(x, proj_id);
% or Create the Spot operator for ASTRA using the GPU.
% A = (dom_width/n)*opTomo('cuda', proj_geom, vol_geom);

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
