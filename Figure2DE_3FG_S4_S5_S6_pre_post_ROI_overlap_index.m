% Determine ROI overlap index between pre- and post-window intrinsic signal imaging after
% finding the best translation.
% Calculate shift_x and shift_y as the shift of center of mass for each area.
%
% Analysis code for Figures 2D, 2E, 3F, 3G, S4, S5, and S6
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.


parentdir = 'YOUR_DATA_DIRECTORY';

areas = {'A1','AAF','VAF','A2'};
colors = {[1 0 0],[0 0.5 0], [0.6 0 0.6],[0 0 1]};

% dilate masks to allow some variability
dilate_masks = 1;
dilate_width = 50; % micrometers.
dilate_width_pix = round(dilate_width/(2300/717));

%% mice

mouse = []; date_pre = []; date_post = []; window_DV = []; window_AP = []; n = 1;

mouse{n} = 'MOUSE1'; date_pre{n} = 'DATE1_pre'; date_post{n} = 'DATE1_post'; window_DV{n} = 1.85; window_AP{n} = 2.76; n = n+1; % dimensions are determined for each mouse
mouse{n} = 'MOUSE2'; date_pre{n} = 'DATE2_pre'; date_post{n} = 'DATE2_post'; window_DV{n} = 2.19; window_AP{n} = 2.93; n = n+1;
% add as many mice as you have



%% go through mice

overlapInd = nan(length(mouse),length(areas)); 
shift_x = nan(length(mouse),length(areas)); % P->A
shift_y = nan(length(mouse),length(areas)); % D->V
for m = 1:length(mouse)

    dir_pre = dir([parentdir filesep date_pre{m} filesep mouse{m} '_*']);
    dir_pre = dir_pre([dir_pre.isdir] == 1);
    if length(dir_pre)~=1
        warning(['multiple folders for ' mouse{m} '. Use the first folder']);
        dir_pre = dir_pre(1);
    end
    dir_post = dir([parentdir filesep date_post{m} filesep mouse{m} '_*']);
    dir_post = dir_post([dir_post.isdir] == 1);
    if length(dir_post)~=1
        warning(['multiple folders for ' mouse{m} '. Use the first folder']);
        dir_post = dir_post(1);
    end

    % load ROIs
    load([dir_pre.folder filesep dir_pre.name filesep 'combined_area_borders.roi4'], '-mat');
    areamasks_pre = mask_combined_area;
    areamasklabels_pre = mask_combined_area_labels;

    load([dir_post.folder filesep dir_post.name filesep 'combined_area_borders.roi4'], '-mat');
    areamasks_post = mask_combined_area;
    areamasklabels_post = mask_combined_area_labels;
    % note that A->P corresponds to rows (dim1), and D->V corresponds to columns (dim2) in masks.

    % window border
    load([dir_post.folder filesep dir_post.name filesep 'window_edge_ROI.roi2'], '-mat'); 
    window_mask = polygon.ROI_mask{1};


    % find the best translation that results in the largest overlaps between area masks.
    % This corresponds to the best translation in 2-dim crosscorrelation. 

    for a = 1:length(areas)
        if ~isempty(areamasks_post{a})
            width = size(areamasks_post{a},2);
            height = size(areamasks_post{a},1);
            break
        end
    end

    CorrCoef_sum = zeros(2*height-1, 2*width-1);
    for a = 1:length(areas)
        currmask_pre = areamasks_pre{strmatch(areas{a},areamasklabels_pre)};
        currmask_post = areamasks_post{strmatch(areas{a},areamasklabels_post)};

        if dilate_masks == 1
            se = strel('disk', dilate_width_pix);
            currmask_pre = imdilate(currmask_pre,se);
            currmask_post = imdilate(currmask_post,se);
        end

        CorrCoef_sum = CorrCoef_sum + xcorr2(double(currmask_pre),double(currmask_post));
    end

    % Find global maximum overlap across all pairs
    [~, idx] = max(CorrCoef_sum(:));
    [peakY, peakX] = ind2sub(size(CorrCoef_sum), idx);

    % Convert peak location to translation of post relative to pre
    dy_prepost = peakY - height;   % row shift (down positive)
    dx_prepost = peakX - width;   % column shift (right positive)

    tform = transltform2d(dx_prepost, dy_prepost); % dx: D->V; dy: A->P; theta: degree, clockwise

    for a = 1:length(areas)
        currmask_pre = areamasks_pre{strmatch(areas{a},areamasklabels_pre)};
        currmask_post = areamasks_post{strmatch(areas{a},areamasklabels_post)};
        
        [~, R_post] = imwarp(currmask_post,tform);
        R_pre = imref2d(size(currmask_pre));

        xLimits = [min(R_pre.XWorldLimits(1), R_post.XWorldLimits(1)), ...
            max(R_pre.XWorldLimits(2), R_post.XWorldLimits(2))];
        yLimits = [min(R_pre.YWorldLimits(1), R_post.YWorldLimits(1)), ...
            max(R_pre.YWorldLimits(2), R_post.YWorldLimits(2))];

        Rcombined = imref2d(round([diff(yLimits), diff(xLimits)]), xLimits, yLimits);

        % Warp pre image into combined coordinates.
        tform_identical = transltform2d(0, 0);
        currmask_pre_combined = imwarp(currmask_pre, tform_identical, 'OutputView', Rcombined);
        currmask_post_combined = imwarp(currmask_post, tform, 'OutputView', Rcombined, 'InterpolationMethod', 'nearest');

        if dilate_masks == 1
            se = strel('disk', dilate_width_pix);
            currmask_pre_combined = imdilate(currmask_pre_combined,se);
            currmask_post_combined = imdilate(currmask_post_combined,se);
        end

        % warp window border 
        windowmask_combined = imwarp(window_mask, tform, 'OutputView', Rcombined, 'InterpolationMethod', 'nearest');

        % remove masks outside the window border
        currmask_pre_combined(~windowmask_combined) = 0;
        currmask_post_combined(~windowmask_combined) = 0;

        % calculate overlap index as Dice index
        sum_masks = currmask_pre_combined+currmask_post_combined;
        overlapInd(m,a) = 2*numel(find(sum_masks==2))/(numel(find(currmask_pre_combined==1))+numel(find(currmask_post_combined==1)));

        % % alternatively, use Jaccard index?
        % overlapInd(m,a) = numel(find(sum_masks==2))/(numel(find(currmask_pre_combined==1))+numel(find(currmask_post_combined==1))-numel(find(sum_masks==2)));
        
        % calculate the shift in the center of mass of the area.
        % note that the definition of x and y changes here to P->A and D->V
        stats_pre = regionprops(currmask_pre_combined,'Centroid');
        stats_post = regionprops(currmask_post_combined,'Centroid');
        shift_x(m,a) = -(stats_post.Centroid(2)-stats_pre.Centroid(2));
        shift_y(m,a) = stats_post.Centroid(1)-stats_pre.Centroid(1);

    end
end


%% quantification

window_DV_all = cell2mat(window_DV);
window_AP_all = cell2mat(window_AP);
overlapInd_all = overlapInd;
windowarea_all = pi.*(window_DV_all/2).*(window_AP_all/2);
shiftdist_all = sqrt(shift_x.^2+shift_y.^2)*(2300/717);


h1 = figure('position',[400,100,340,300]);
hold on;
for a = 1:length(areas)
    scatter(window_DV_all,overlapInd_all(:,a),'markeredgecolor',colors{a});
    pl = polyfit(window_DV_all(~isnan(window_DV_all)),overlapInd_all(~isnan(window_DV_all),a),1);
    line([1.4 2.8],[1.4 2.8]*pl(1)+pl(2),'color',colors{a});
end
xlabel('Window D-V size (mm)');
ylabel('Overlap Ind');
ylim([0 1]);
xlim([1.2 3]);

h2 = figure('position',[700,100,400,300]);
hold on;
for a = 1:length(areas)
    scatter(window_AP_all,overlapInd_all(:,a),'markeredgecolor',colors{a});
    pl = polyfit(window_AP_all(~isnan(window_AP_all)),overlapInd_all(~isnan(window_AP_all),a),1);
    line([2.3 3.7],[2.3 3.7]*pl(1)+pl(2),'color',colors{a});
end
xlabel('Window A-P size (mm)');
ylabel('Overlap Ind');
ylim([0 1]);

h3 = figure('position',[1000,100,400,300]);
hold on;
for a = 1:length(areas)
    scatter(windowarea_all,overlapInd_all(:,a),'markeredgecolor',colors{a});
    pl = polyfit(windowarea_all(~isnan(windowarea_all)),overlapInd_all(~isnan(windowarea_all),a),1);
    line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color',colors{a});
end
xlabel('Window area (mm^2)');
ylabel('Overlap Ind');
ylim([0 1]);
xlim([2.5 7.5]);

[R_OI_A1 p_OI_A1] = corrcoef(windowarea_all,overlapInd_all(:,1));
[R_OI_AAF p_OI_AAF] = corrcoef(windowarea_all,overlapInd_all(:,2));
[R_OI_VAF p_OI_VAF] = corrcoef(windowarea_all,overlapInd_all(:,3));
[R_OI_A2 p_OI_A2] = corrcoef(windowarea_all,overlapInd_all(:,4));


%% alternatively, give one value per mouse


h1 = figure('position',[400,100,340,300]);
hold on;
scatter(window_DV_all,mean(overlapInd_all,2),'markeredgecolor','k');
pl = polyfit(window_DV_all(~isnan(window_DV_all)),mean(overlapInd_all(~isnan(window_DV_all),:),2),1);
[R_DV p_DV] = corrcoef(window_DV_all(~isnan(window_DV_all)),mean(overlapInd_all(~isnan(window_DV_all),:),2));
line([1.4 2.8],[1.4 2.8]*pl(1)+pl(2),'color','r');
xlabel('Window D-V size (mm)');
ylabel('Overall Overlap Ind');
ylim([0 1]);
xlim([1.2 3]);

h2 = figure('position',[700,100,340,300]);
hold on;
scatter(window_AP_all,mean(overlapInd_all,2),'markeredgecolor','k');
pl = polyfit(window_AP_all(~isnan(window_AP_all)),mean(overlapInd_all(~isnan(window_AP_all),:),2),1);
[R_AP p_AP] = corrcoef(window_AP_all(~isnan(window_AP_all)),mean(overlapInd_all(~isnan(window_AP_all),:),2));
line([2 3.8],[2 3.8]*pl(1)+pl(2),'color','r');
xlabel('Window A-P size (mm)');
ylabel('Overall Overlap Ind');
ylim([0 1]);
xlim([2 3.8]);


% instead, take a one "dominating" window size.
% single explaining variables like max([window_DV_all, window_AP_all*0.8]) may work better.
a = window_DV_all(:); b = window_AP_all(:); X = mean(overlapInd_all,2);

% grid search for the best 'explaining variable'

r_vals = 0:0.01:2;
R2 = nan(size(r_vals));

for i = 1:length(r_vals)
    z = max(a, r_vals(i)*b);
    if std(z) ~= 0
        r = corr(z, X, 'Rows','complete');
        R2(i) = r.^2;
    end
end

[~, i_best] = max(R2);
r_best = r_vals(i_best);

effective_size = max(a, r_best*b);

h3 = figure('position',[100 100 250 250]);
plot(r_vals,R2(:),'-k');
xlim([0.5 1]);
set(gca,'box','off');
xlabel('AP scaling factor');
ylabel('R2');

h4 = figure('position',[1100,100,340,300]);
hold on;
scatter(effective_size,mean(overlapInd_all,2),'markeredgecolor','k');
pl = polyfit(effective_size,mean(overlapInd_all,2),1);
[R p] = corrcoef(effective_size,mean(overlapInd_all,2));
line([1.2 3],[1.2 3]*pl(1)+pl(2),'color','r');
xlabel('Best explaining variable (mm)');
ylabel('Overall Overlap Ind');
xlim([1.2 3]);
ylim([0 1]);

% or use window area
h5 = figure('position',[1400,100,340,300]);
hold on;
scatter(windowarea_all,mean(overlapInd_all,2),'markeredgecolor','k');
pl = polyfit(windowarea_all,mean(overlapInd_all,2),1);
[R p] = corrcoef(windowarea_all,mean(overlapInd_all,2));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
xlabel('Window area (mm2)');
ylabel('Overall Overlap Ind');
ylim([0 1]);
xlim([2.5 7.5]);


%% quantification

mean_overlapInd = mean(overlapInd_all,2);
h = figure('position',[100 400 200 250]); hold on;
plotSpread_modified(mean_overlapInd,'distributionMarkers','o','distributionColors','k');
line([0.5 1.5],mean(mean_overlapInd)*[1 1],'color','r','linewidth',1.5);
ylim([0 1]);
set(gca,'xtick',[]);
ylabel('Overlap Index');


% separately for individual areas
hh = figure('position',[400 400 400 250]);
hold on;
data = {};
for a = 1:length(areas)
    data{1,a} = overlapInd_all(:,a);
end

plotSpread_modified(data,'distributionMarkers','o','distributionColors',colors);
for a = 1:length(areas)
    line([a-0.4 a+0.4],mean(overlapInd_all(:,a))*[1 1],'color','k','linewidth',1.5);
end

set(gca,'xtick',1:length(areas),'xticklabel',areas);
ylim([0 1]);
ylabel('Overlap Index');


%% plot shift_x and shift_y

h_shift = figure('position',[300 300 300 300]); hold on;

for a = 1:length(areas)
    % scatter
    scatter(shift_x(:,a)*(2300/717),shift_y(:,a)*(2300/717),'markerfacecolor','none','markeredgecolor',colors{a},'markeredgealpha',0.3); % convert to um

end

xlabel('X shift (um; P->A)');
ylabel('Y shift (um; D->V)');
set(gca,'ydir','reverse');
xlim([-500 500]);
ylim([-500 500]);
line([-500 500],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
line([0 0],[-500 500],'linestyle',':','color',[0.5 0.5 0.5]);

for a = 1:length(areas)
    plot(mean(shift_x(:,a)*(2300/717)),mean(shift_y(:,a)*(2300/717)),'marker','x','markerfacecolor','none','markeredgecolor',([0 0 0] + colors{a})*2/3,'markersize',12); % make markers darker
end


% x and y axes separately
h_shift_x = figure('position',[600 300 300 200]); hold on;
for a = 1:length(areas)
    scatter(a*ones(size(shift_x,1),1)+0.1*randn(size(shift_x,1),1),shift_x(:,a)*(2300/717),'markerfacecolor','none','markeredgecolor',colors{a},'markeredgealpha',0.5);
    line([a-0.3 a+0.3],mean(shift_x(:,a)*(2300/717))*[1 1],'color','k','linewidth',1.5);

end
line([0 5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
ylim([-500 500]);
set(gca,'xtick',1:4,'xticklabel',areas);
ylabel('X shift (um; P->A)');

h_shift_y = figure('position',[900 300 300 200]); hold on;
for a = 1:length(areas)
    scatter(a*ones(size(shift_y,1),1)+0.1*randn(size(shift_y,1),1),shift_y(:,a)*(2300/717),'markerfacecolor','none','markeredgecolor',colors{a},'markeredgealpha',0.5);
    line([a-0.35 a+0.35],mean(shift_y(:,a)*(2300/717))*[1 1],'color','k','linewidth',1.5);
end
line([0 5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
ylim([-500 500]);
set(gca,'xtick',1:4,'xticklabel',areas,'ydir','reverse');
ylabel('Y shift (um; D->V)');

p_A1_x = signrank(shift_x(:,strcmp(areas,'A1')))*4;
p_AAF_x = signrank(shift_x(:,strcmp(areas,'AAF')))*4;
p_VAF_x = signrank(shift_x(:,strcmp(areas,'VAF')))*4;
p_A2_x = signrank(shift_x(:,strcmp(areas,'A2')))*4;
p_A1_y = signrank(shift_y(:,strcmp(areas,'A1')))*4;
p_AAF_y = signrank(shift_y(:,strcmp(areas,'AAF')))*4;
p_VAF_y = signrank(shift_y(:,strcmp(areas,'VAF')))*4;
p_A2_y = signrank(shift_y(:,strcmp(areas,'A2')))*4;



% plot absolute shift distance against window size

h1 = figure('position',[400,100,340,300]);
hold on;
scatter(window_DV_all,mean(shiftdist_all,2),'markeredgecolor','k');
pl = polyfit(window_DV_all(~isnan(window_DV_all)),mean(shiftdist_all(~isnan(window_DV_all),:),2),1);
[R_DV p_DV] = corrcoef(window_DV_all(~isnan(window_DV_all)),mean(shiftdist_all(~isnan(window_DV_all),:),2));
line([1.4 2.8],[1.4 2.8]*pl(1)+pl(2),'color','r');
xlabel('Window D-V size (mm)');
ylabel('Shift distance (um)');
ylim([0 400]);
set(gca,'ytick',0:100:400);
xlim([1.2 3]);

h2 = figure('position',[700,100,340,300]);
hold on;
scatter(window_AP_all,mean(shiftdist_all,2),'markeredgecolor','k');
pl = polyfit(window_AP_all(~isnan(window_AP_all)),mean(shiftdist_all(~isnan(window_AP_all),:),2),1);
[R_AP p_AP] = corrcoef(window_AP_all(~isnan(window_AP_all)),mean(shiftdist_all(~isnan(window_AP_all),:),2));
line([2 3.8],[2 3.8]*pl(1)+pl(2),'color','r');
xlabel('Window A-P size (mm)');
ylabel('Shift distance (um)');
ylim([0 400]);
xlim([2 3.8]);
set(gca,'ytick',0:100:400,'xtick',2.5:0.5:3.5);


% instead, take a one "dominating" window size.
% single explaining variables like max([window_DV_all, window_AP_all*0.8]) may work better.
a = window_DV_all(:); b = window_AP_all(:); X = mean(shiftdist_all,2);

% grid search for the best 'explaining variable'

r_vals = 0:0.01:2;
R2 = nan(size(r_vals));

for i = 1:length(r_vals)
    z = max(a, r_vals(i)*b);
    if std(z) ~= 0
        r = corr(z, X, 'Rows','complete');
        R2(i) = r.^2;
    end
end

[~, i_best] = max(R2);
r_best = r_vals(i_best);

effective_size = max(a, r_best*b);

h3 = figure('position',[100 100 250 250]);
plot(r_vals,R2(:),'-k');
xlim([0.5 1]);
set(gca,'box','off');
xlabel('AP scaling factor');
ylabel('R2');

h4 = figure('position',[1100,100,340,300]);
hold on;
scatter(effective_size,mean(shiftdist_all,2),'markeredgecolor','k');
pl = polyfit(effective_size,mean(shiftdist_all,2),1);
[R p] = corrcoef(effective_size,mean(shiftdist_all,2));
line([1.4 3.2],[1.4 3.2]*pl(1)+pl(2),'color','r');
xlabel('Best explaining variable (mm)');
ylabel('Shift distance (um)');
xlim([1.2 3]);
ylim([0 400]);
set(gca,'ytick',0:100:400);

% or use window area
h5 = figure('position',[1400,100,340,300]);
hold on;
scatter(windowarea_all,mean(shiftdist_all,2),'markeredgecolor','k');
pl = polyfit(windowarea_all,mean(shiftdist_all,2),1);
[R p] = corrcoef(windowarea_all,mean(shiftdist_all,2));
line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color','r');
xlabel('Window area (mm2)');
ylabel('Shift distance (um)');
xlim([2.5 7.5]);
ylim([0 400]);
set(gca,'ytick',0:100:400);


% plot individual areas
h6 = figure('position',[1700,100,400,300]);
hold on;
for a = 1:length(areas)
    scatter(windowarea_all,shiftdist_all(:,a),'markeredgecolor',colors{a});
    pl = polyfit(windowarea_all(~isnan(windowarea_all)),shiftdist_all(~isnan(windowarea_all),a),1);
    line([2.5 7.5],[2.5 7.5]*pl(1)+pl(2),'color',colors{a});
end
xlabel('Window area (mm^2)');
ylabel('Shift distance (um)');
ylim([0 450]);
xlim([2.5 7.5]);
set(gca,'ytick',0:100:400);

[R_shift_A1 p_shift_A1] = corrcoef(windowarea_all,shiftdist_all(:,1));
[R_shift_AAF p_shift_AAF] = corrcoef(windowarea_all,shiftdist_all(:,2));
[R_shift_VAF p_shift_VAF] = corrcoef(windowarea_all,shiftdist_all(:,3));
[R_shift_A2 p_shift_A2] = corrcoef(windowarea_all,shiftdist_all(:,4));
