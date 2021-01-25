
%% Denoising Iterative Hessian Sketching (Denoising-IHS)
%| [Ref.] DOI: 
%| Copyright 25-01-2021, Alessandro Perelli, Technical University of Denmark (DTU)

%% ******************  Description  **********************
%| This script upload the spectral CT data, post-process the data 
%| to insert gaps and generate the ASTRA operator.
%| As option, the SIRT reconstruction can be performed.
%| option_recon = 1;

%| This code requires the ASTRA Toolbox installed:
%| https://github.com/astra-toolbox/astra-toolbox
%% ********************************************************

%% ******************  Requirements  **********************
%| the dataset can be downloaded at the following link:

%| W. van Aarle, W. J. Palenstijn, J. Cant, E. Janssens, F. Bleichrodt, 
%| A. Dabravolski, J. De Beenhouwer, K. J. Batenburg, and J. Sijbers, 
%| “Fast and Flexible X-ray Tomography Using the ASTRA Toolbox”, 
%| Optics Express, 24(22), 25129-25147, (2016), 
%| http://dx.doi.org/10.1364/OE.24.025129

%% ****************************************************

clear, close all;

opt_recon = 0;

% This script is assuming single energy 3D volumetric data. 
a  = h5info('HER_acq19-02-13-16-00-40.h5');
E  = 128;
nr = 640;
nc = 100;
nc_p = 3700;

data      = zeros(nr,nc_p,E);
data_flat = zeros(nr,nc,E);

for i = 1:E
    data_flat(:,:,i) = h5read('HER_acq19-02-13-16-00-40.h5',['/' a.Datasets(i).Name]);
    data(:,:,i)=  h5read('HER_acq19-02-13-16-01-00.h5',['/' a.Datasets(i).Name]);
end

flat = mean(data_flat,2);
for i=1:360
    data_n(:,i,:) = mean(data(:,(5+(i-1)*10+1):(5+(i)*10+1),:),2);
end

for Ech = 1:E
    data_E = squeeze(data_n(:,:,Ech));
    flat_E = squeeze(flat(:,:,Ech));
    data_E = [flat_E data_E]';
    % fill 4 pixels gap between 2 multix modules by neigbour interpolation
    data_E = insert_gap(data_E,128*4,0,4);
    data_E = insert_gap(data_E,128*3,0,4);
    data_E = insert_gap(data_E,128*2,0,4);
    data_E = insert_gap(data_E,128,0,4);

    data_E = data_E';
    
    for i=1:size(data_E,2)-1
        %     tomo(:,:,i) = att_pos(data_E(:,:,i+1),data_E(:,:,1),0,0);
        tomo(:,i) = -log(data_E(:,i+1)./data_E(:,1));
    end
    %     Slicer(tomo) %see data slice by slice
    % tomo = permute(tomo,[1 3 2]); %reshaping to fit astra format
    % tomo(:,:,575:end)=[]; % cutting sample holder (might need to be reset)

    %% Reconstruction is with fan beams geometry (translating line 1D-detector). 
    %  For multi energy please loop it to fit the data size (or optimize it differently)
    if opt_recon == 1 

        %setting up astra geometry
        det_spacing_x = 0.08025; %pixel size x
        det_col_count = size(tomo,1); %pixel number y
        source_origin = 81.66; %source to axis of rotation
        origin_det = 114.19 - source_origin; %source to axis of rotation

        %rescaling to make astra work (if you find why this is necessary or a better way to do it, please let
        %me know)
        M = (origin_det+source_origin)/source_origin;
        scale_x = M/det_spacing_x;
        % scale_y = M/det_spacing_y;
        source_origin = scale_x*source_origin;
        origin_det = scale_x*origin_det;
        det_spacing_x = det_spacing_x*scale_x;
        % det_spacing_y = det_spacing_y*scale_y;

        % vol_geom = astra_create_vol_geom(size(tomo,1),size(tomo,1),size(tomo,3));
        vol_geom = astra_create_vol_geom(size(tomo,1),size(tomo,1));

        src_shift = det_spacing_x*0; % source shift
        det_shift = 4.41*det_spacing_x; % detector shift
        det_tilt = 0; % detector tilt (radians)
        center_shift = 0; % shift of axis of rotation
        tomo_sh = circshift(tomo,center_shift,1);
        tomofl = flipud(tomo_sh); % flipping sinogram 
        angles = linspace(0,2*pi,size(tomo_sh,2)+1); angles = angles(1:size(tomo_sh,2)); 
        vectors = fan_vec_custom(angles,source_origin,origin_det,det_spacing_x,src_shift,det_shift,det_tilt);
        
        proj_geom = astra_create_proj_geom('fanflat_vec', det_col_count, vectors);
        
        attenrad = tomo_sh';
        sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, attenrad);

        disp('Performing Reconstruction'), tic
        k = 50;
        rec_id = astra_mex_data2d('create', '-vol', vol_geom);

        % cfg = astra_struct('CGLS3D_CUDA');
        cfg = astra_struct('SIRT');
        % cfg = astra_struct('FDK_CUDA');
        cfg.option.MinConstraint = 0;
        cfg.ReconstructionDataId = rec_id;
        cfg.ProjectionDataId = sinogram_id;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        proj_id = astra_create_projector('linear', proj_geom, vol_geom);
        cfg.ProjectorId      = proj_id;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        alg_id = astra_mex_algorithm('create', cfg);
        astra_mex_algorithm('iterate', alg_id, k);
        % astra_mex_algorithm('run', alg_id);
        reconstr(:,:,Ech) = astra_mex_data2d('get', rec_id); toc,
        astra_mex_data2d('delete', rec_id);
        astra_mex_data2d('delete', sinogram_id);
        astra_mex_data2d('delete', alg_id);
        astra_mex_data2d('clear');
    end
end

if opt_recon == 1
    reconstr(reconstr<0)=0;
    % Slicer(reconstr)
    % volumeViewer(reconstr)
    Eb = linspace(20,160,128);
    figure,
    plot(Eb,squeeze(reconstr(300,325,:))), hold on
    plot(Eb,squeeze(reconstr(325,357,:))), hold on
    plot(Eb,squeeze(reconstr(355,331,:))), hold on
    plot(Eb,squeeze(reconstr(331,297,:))), hold on
end

save('tomo.mat', 'tomo', 'data_E');
