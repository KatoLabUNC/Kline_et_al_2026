% Show  heatmap of harmonics responses averaged across animals,
% centered and aligned around the auditory cortex center.
%
% Analysis code for Figures 1D and 4C
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.


scalefactor = 2.3/717; % mm/px

parentdir = 'YOUR_DATA_DIRECTORY';

areas = {'VAF','AAF','A2','A1'}; % changed the order so A2 and A1 come to the top layers of the plot

% align_center = 'AC_center'; % center of the entire Aud

deblurring = 1;
filter_size = 10;

display_method = 'Zscore_from_mu';
% display_method = 'raw_percent';


%% mice
n = 1;
mouse = []; 
date = [];
expdir = [];
rotation = []; % radian. determined based on the ridge. (based on Narayanan 2023 data, AP axis is roughly parallel to the ridge around AC.)

% all mice
mouse{n} = 'MOUSE1'; date{n} = 'DATE1'; expdir{n} = 'EXPDIR1'; rotation{n} = 0.232; n = n+1; % rotation values are determined for each mouse
mouse{n} = 'MOUSE2'; date{n} = 'DATE2'; expdir{n} = 'EXPDIR2'; rotation{n} = 0.206; n = n+1;
mouse{n} = 'MOUSE3'; date{n} = 'DATE3'; expdir{n} = 'EXPDIR3'; rotation{n} = 0.339; n = n+1;
% add as many mice as you have



%%

harmonicsMap_all = [];
for m = 1:length(mouse)
    mousedir = [parentdir filesep date{m} filesep expdir{m}];

    % get all AC masks to determine the AC center
    load([mousedir filesep 'combined_area_borders.roi4'],'-mat');
    masks_include = find(ismember(mask_combined_area_labels,{'A1','A2','VAF','AAF'}));
    combined_mask = zeros(size(mask_combined_area{1}));
    for a = 1:length(masks_include)
        combined_mask = combined_mask+mask_combined_area{masks_include(a)};
    end
    AC_center = [(find(sum(combined_mask,2)>0, 1 )+find(sum(combined_mask,2)>0, 1, 'last' ))/2 ...
        (find(sum(combined_mask,1)>0, 1 )+find(sum(combined_mask,1)>0, 1, 'last' ))/2];
    align_center_coor = AC_center; % y,x = (A->P, D->V) from top-left


    A1_ind = find(strcmp('A1',mask_combined_area_labels));
    A1_mask = mask_combined_area{A1_ind};
    A2_ind = find(strcmp('A2',mask_combined_area_labels));
    A2_mask = mask_combined_area{A2_ind};
    VAF_ind = find(strcmp('VAF',mask_combined_area_labels));
    VAF_mask = mask_combined_area{VAF_ind};
    AAF_ind = find(strcmp('AAF',mask_combined_area_labels));
    AAF_mask = mask_combined_area{AAF_ind};

    % apply rotation and scaling
    for a = 1:length(areas)
        eval(['curr_mask = ' areas{a} '_mask;']);

        curr_border = mask2poly(curr_mask,'exact'); % (x,y) in the direction of D->V, A->P,
        curr_border(1,:) = [];

        % convert to millimeter
        curr_border = scalefactor*(curr_border-[align_center_coor(2) align_center_coor(1)]);
        curr_border(:,2) = -curr_border(:,2); % change to P->A

        % now, the AC centroid is at [0 0], so rotation can be applied around it
        % rotation matrix is
        % [cos sin 0;
        % -sin cos 0;
        %   0   0  1]

        A = [cos(-rotation{m}) sin(-rotation{m}) 0; ...
            -sin(-rotation{m}) cos(-rotation{m}) 0; ...
            0         0        1];
        tform = maketform('affine',A);
        curr_border = tformfwd(tform,curr_border);

        curr_border = fliplr(curr_border); % (x,y) = (A->P, D->V)

        eval([areas{a} '_borders_all{m} = curr_border;']);

    end

    % re-center the auditory cortex. 
    % Now, A-P and D-V axes are aligned with x and y.
    % (I could have simply rotate before shifting.)

    borders_all = [];
    for a = 1:length(areas)
       eval(['borders_all = [borders_all; ' areas{a} '_borders_all{m}];']);
    end
    AC_center = [(min(borders_all(:,1))+max(borders_all(:,1)))/2, ...
        (min(borders_all(:,2))+max(borders_all(:,2)))/2]; % in millimeter
    align_center_coor2 = AC_center/scalefactor; % (x,y) = (A->P, D->V) % pixels
    
    
    % load harmonics response image

    % find the avg image corresponding to this sound
    files_curr = dir([mousedir filesep 'harmonics_img_stack__div_%change*.tif']);
    
    % if multiple, take the first one
    tiff_fullfilename = [mousedir filesep files_curr(1).name];

    im = [];
    im(:,:)=double(imread([tiff_fullfilename],'tiff',1));
    if min(im(:))>=2^15
        im = im-2^15;
    end
    
    % deblur image
    if deblurring == 1
        % deblur the image.(Issa et al., Neuron 83, 944)
        % point spread function. assume gaussian with 200um sigma.
        % Currently, 200um corresponds to 68 pixels.
        PSF = fspecial('gaussian', 400,68);
        
        im = -deconvlucy(-im+1,PSF,3)+1;
    end

    % apply filter to smoothen images
    H = fspecial('disk', filter_size);
    FilteredImg = imfilter(im, H, 'replicate');
    
    switch display_method
        case 'Zscore_from_mu'
            % the fit with larger mu should be the baseline
            reset(RandStream.getGlobalStream);
            obj = gmdistribution.fit(im(im~=0), 2);
            [Baseline, ind] = max(obj.mu);
            BaselineSigma = obj.Sigma(1,1,ind);

            im_sound = (FilteredImg-Baseline)/sqrt(BaselineSigma);
        
        case 'raw_percent'
            im_sound = FilteredImg;
    end


    % apply rotation and translation
    outH = 2048; outW = 2048;
    cx = (outW+1)/2;
    cy = (outH+1)/2;

    Rout = imref2d([outH outW]);    % fixed output grid
    Rin  = imref2d(size(im_sound));        % input grid (intrinsic coords)

    % 1) translate (x,y) -> center of output canvas
    tx = cx - align_center_coor(2);
    ty = cy - align_center_coor(1);
    T_translate = [1 0 0;
                   0 1 0;
                   tx ty 1];

    % 2) rotate about output center
    R = [ cos(rotation{m})  sin(rotation{m})  0;
        -sin(rotation{m})  cos(rotation{m})  0;
        0         0         1];

    T_toCenter   = [1 0 0; 0 1 0; cx cy 1];
    T_fromCenter = [1 0 0; 0 1 0;  -cx  -cy 1];
    T_rotateAboutCenter = T_fromCenter * R * T_toCenter;
    
    % 3) translate again so that the center is defined along precise A-P and D-V axes.
    tx2 = align_center_coor2(2);
    ty2 = align_center_coor2(1);
    T_translate2 = [1 0 0;
        0 1 0;
        tx2 ty2 1];    

    % Apply translate first, then rotate about the canvas center
    T_total = T_translate* T_rotateAboutCenter * T_translate2;

    tform = affine2d(T_total);

    im_sound_centered_rotated = imwarp(im_sound, Rin, tform, ...
        'OutputView', Rout, ...
        'Interp', 'linear', ...
        'FillValues', 0);

    harmonicsMap_all(:,:,m) = im_sound_centered_rotated;
end


%% plot

harmonics_avgMap = nanmean(harmonicsMap_all,3);

% rotate by 90 degree
harmonics_avgMap = imrotate(harmonics_avgMap, -90);

h = figure('position',[100 100 300 300]);
imagesc(harmonics_avgMap,[-5 0]);
colormap(flipud(jet));
set(gca,'xtick',[],'ytick',[]);
xlim([-1 1]/scalefactor+1024.5);
ylim([-1 1]/scalefactor+1024.5);
tickpositions = (-1:0.5:1)/scalefactor+1024.5;
set(gca,'xtick',tickpositions,'ytick',tickpositions,'xticklabel',-1:0.5:1, 'yticklabel',-1:0.5:1, 'box','off');
xlabel('Anterior from AC Center (mm)');
ylabel('Ventral from AC Center (mm)');


switch display_method
    case 'Zscore_from_mu'
        % optionally, adjust lookup table
        climits = [-5 0];
        set(findobj(h,'type','axes'),'clim',climits);

        hh = figure;
        %set the color map for black and all colors
        jet2 = jet(100);
        A = 1:100;
        colorbar = repmat(flipud(A'),1,10);
        imshow(colorbar, jet2);
        title(['-zscore min:' num2str(-climits(2)) ' max:' num2str(-climits(1))]);
    case 'raw_percent'
        % optionally, adjust lookup table
        climits = [-0.15 0];
        set(findobj(h,'type','axes'),'clim',climits);

        hh = figure;
        %set the color map for black and all colors
        jet2 = jet(100);
        A = 1:100;
        colorbar = repmat(flipud(A'),1,10);
        imshow(colorbar, jet2);
        title(['percent min:' num2str(-climits(2)) ' max:' num2str(-climits(1))]);
end


