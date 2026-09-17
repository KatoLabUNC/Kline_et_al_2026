% Plot relative coordinates of auditory cortical areas in intrinsic imaging mapping data.
% Plot A1, A2, AAF, and VAF relative to the AC Center.
%
% Analysis code for Figures 1D, 3B-E, and 4C
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.


scalefactor = 2.3/717; % scaling factor. mm per pixels

areas = {'VAF','AAF','A2','A1'}; % changed the order so A2 and A1 come to the top layers of the plot
area_colors = {[0.6 0 0.6],[0 0.5 0],[0 0 1],[1 0 0]};

parentdir = 'YOUR_DATA_DIRECTORY';



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

A1_borders_all = [];
A2_borders_all = [];
VAF_borders_all = [];
AAF_borders_all = [];
for m = 1:length(mouse)
    mousedir = [parentdir filesep date{m} filesep expdir{m}];


    % load all area masks
    load([mousedir filesep 'combined_area_borders.roi4'],'-mat');
    masks_include = find(ismember(mask_combined_area_labels,{'A1','A2','VAF','AAF'}));
    combined_mask = zeros(size(mask_combined_area{1}));
    for k = 1:numel(masks_include)
        idx = masks_include(k);
        combined_mask = combined_mask + mask_combined_area{idx};
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

    borders_all = [];
    for a = 1:length(areas)
       eval(['borders_all = [borders_all; ' areas{a} '_borders_all{m}];']);
    end
    AC_center = [(min(borders_all(:,1))+max(borders_all(:,1)))/2, ...
        (min(borders_all(:,2))+max(borders_all(:,2)))/2];
    align_center_coor2 = AC_center; % (x,y) = (A->P, D->V) % in millimeter
    
    for a = 1:length(areas)
        eval([areas{a} '_borders_all{m} = ' areas{a} '_borders_all{m} - align_center_coor2;']);
    end

end


%% plot polygons overlaid

% border coordinates are (x,y) in the direction of (A->P, D->V), with the AC center as zero
h = figure('position',[100 100 300 300]);
hold on;
set(gca,'Ydir','reverse','xtick',-1:0.5:1,'ytick',-1:0.5:1);

for m = 1:length(mouse)
    for a = 1:length(areas)

        eval(['curr_border = ' areas{a} '_borders_all{m};']);

        if isempty(curr_border)
            continue
        end

        plot(curr_border(:,1),curr_border(:,2),'color',area_colors{a});
    end
end

axis equal
xlim([-1 1]);
ylim([-1 1]);
line([-1 1],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
line([0 0],[-1 1],'linestyle',':','color',[0.5 0.5 0.5]);

xlabel('Anterior from AC center (mm)');
ylabel('Ventral from AC center (mm)');



%% plot probability distribution

% change origin to (x,y) = [-1.5, -1] (dorso-posterior corner)
counts = zeros(1,length(areas)); % mouse count
for m = 1:length(mouse)
    for a = 1:length(areas)
        if ~isempty(eval([areas{a} '_borders_all{m}']))
            eval([areas{a} '_borders_all2{m} = (' areas{a} '_borders_all{m} + [1.5 1])*1000;']); % shift origin and convert to um
            counts(a) = counts(a)+1;
        else
            eval([areas{a} '_borders_all2{m} = [];']);
        end
    end
end


% convert to mask
A1_mask_prob = zeros(5000,2500); %(P->A,D->V)
VAF_mask_prob = zeros(5000,2500);
AAF_mask_prob = zeros(5000,2500);
A2_mask_prob = zeros(5000,2500);
for m = 1:length(mouse)
    for a = 1:length(areas)
        try
            eval([areas{a} '_mask_curr = poly2mask(' areas{a} '_borders_all2{m}(:,2),' areas{a} '_borders_all2{m}(:,1),5000,2500);']);
            eval([areas{a} '_mask_prob =' areas{a} '_mask_prob + ' areas{a} '_mask_curr/counts(a);']);
        catch
            continue
        end
    end
end


downsample_factor = 5;

% smoothen
smooth_image = 1;
if smooth_image == 1
    filter_size = 150;
    H = fspecial('disk', filter_size);
end

A1_mask_prob = imfilter(A1_mask_prob, H, 'replicate');
A1_mask_prob = imresize(A1_mask_prob, 1/downsample_factor);
A1_mask_prob(A1_mask_prob<0) = 0;
A2_mask_prob = imfilter(A2_mask_prob, H, 'replicate');
A2_mask_prob = imresize(A2_mask_prob, 1/downsample_factor);
A2_mask_prob(A2_mask_prob<0) = 0;
AAF_mask_prob = imfilter(AAF_mask_prob, H, 'replicate');
AAF_mask_prob = imresize(AAF_mask_prob, 1/downsample_factor);
AAF_mask_prob(AAF_mask_prob<0) = 0;
VAF_mask_prob = imfilter(VAF_mask_prob, H, 'replicate');
VAF_mask_prob = imresize(VAF_mask_prob, 1/downsample_factor);
VAF_mask_prob(VAF_mask_prob<0) = 0;

hh = figure('position',[500 100 300 300]);
hold on;
X = 0:downsample_factor:5000; X(end) = [];
Y = 0:downsample_factor:2500; Y(end) = [];

X = (X-1500)/1000;
Y = (Y-1000)/1000;

% plot contour as line object 
C_A1 = contourc(X,Y,A1_mask_prob',0.199:0.1:1);
k = 1;
while k<size(C_A1,2)
    level = C_A1(1,k);
    npts = C_A1(2,k);
    pts = C_A1(:,k+1:k+npts);
    plot(pts(1,:),pts(2,:), 'color', area_colors{4}); % This may look strange, but contourc output has "header" columns at the beginning of each level, with everything concatenated as a matrix.
    k = k+npts+1;
end
C_VAF = contourc(X,Y,VAF_mask_prob',0.199:0.1:1);
k = 1;
while k<size(C_VAF,2)
    level = C_VAF(1,k);
    npts = C_VAF(2,k);
    pts = C_VAF(:,k+1:k+npts);
    plot(pts(1,:),pts(2,:), 'color', area_colors{1}); 
    k = k+npts+1;
end
C_A2 = contourc(X,Y,A2_mask_prob',0.199:0.1:1);
k = 1;
while k<size(C_A2,2)
    level = C_A2(1,k);
    npts = C_A2(2,k);
    pts = C_A2(:,k+1:k+npts);
    plot(pts(1,:),pts(2,:), 'color', area_colors{3}); 
    k = k+npts+1;
end
C_AAF = contourc(X,Y,AAF_mask_prob',0.199:0.1:1);
k = 1;
while k<size(C_AAF,2)
    level = C_AAF(1,k);
    npts = C_AAF(2,k);
    pts = C_AAF(:,k+1:k+npts);
    plot(pts(1,:),pts(2,:), 'color', area_colors{2}); 
    k = k+npts+1;
end

set(gca,'Ydir','reverse','xtick',-1:3,'ytick',-1:0.5:1.5);
axis equal
xlim([-1 1]);
ylim([-1 1]);
line([-1 1],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
line([0 0],[-1 1],'linestyle',':','color',[0.5 0.5 0.5]);

xlabel('Anterior from AC center (mm)');
ylabel('Ventral from AC center (mm)');


% Only 50% contour
show_contour = 0.5; 

hhh = figure('position',[500 100 300 300]);
hold on;
X = 0:downsample_factor:5000; X(end) = [];
Y = 0:downsample_factor:2500; Y(end) = [];

X = (X-1500)/1000;
Y = (Y-1000)/1000;

% plot contour as line object
C_A1 = contourc(X,Y,A1_mask_prob',show_contour*[1 1]);
C_VAF = contourc(X,Y,VAF_mask_prob',show_contour*[1 1]);
C_A2 = contourc(X,Y,A2_mask_prob',show_contour*[1 1]);
C_AAF = contourc(X,Y,AAF_mask_prob',show_contour*[1 1]);

% the first column is a header
plot(C_A1(1,2:end),C_A1(2,2:end),'color',area_colors{4});
plot(C_VAF(1,2:end),C_VAF(2,2:end),'color',area_colors{1});
plot(C_A2(1,2:end),C_A2(2,2:end),'color',area_colors{3});
plot(C_AAF(1,2:end),C_AAF(2,2:end),'color',area_colors{2});


set(gca,'Ydir','reverse','xtick',-1:3,'ytick',-1:0.5:1.5);
axis equal
xlim([-1 1]);
ylim([-1 1]);
line([-1 1],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
line([0 0],[-1 1],'linestyle',':','color',[0.5 0.5 0.5]);

xlabel('Anterior from AC center (mm)');
ylabel('Ventral from AC center (mm)');


