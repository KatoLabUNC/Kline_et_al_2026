% calculate coincidence-preference index and compare with window size
%
% Analysis code for Figure 5
%
% Mouse identifiers, directory names, and per-mouse parameters shown below
% are placeholders/examples and should be replaced with dataset-specific values.


force_nonresponsive_to_zero = 0; % fix to 0.
include_negative = 1; % fix to 1.

area_threshold = 1; % Include only ROIs with at least this number of responses.
F_threshold = 0; % Only include ROIs with at least this baseline fluorescence.

include_dTind = [2 3 5 6];


% % select the cells that you want to show
% % (choose one of the options)

% showROI = 'All'; % show all the ROIs
% showROI = 'All_accept'; % show all the ROIs which are not rejected
showROI = 'nolabel'; % show all the ROIs without label (no GAD, not rejected)
% showROI = 'GAD'; % show only the GAD cells


parentdir = 'YOUR_DATA_DIRECTORY';


%% enter mouse information/define paths for harmonics_F0 experiment

mouse = []; date = []; layer = []; window_DV = []; window_AP = []; n = 1; 

% %%%%%%%%%%%%%%%%%%%%%%%%%
% % A1 L2/3 pyramidal cells
% %%%%%%%%%%%%%%%%%%%%%%%%%

mouse{n} = 'MOUSE1'; date{n} = 'DATE1'; layer{n}{1} = 'REGIONNAME1'; window_DV{n} = 1.61; window_AP{n} = 2.71;  n = n+1; % window dimensions are determined for each mouse
mouse{n} = 'MOUSE2'; date{n} = 'DATE2'; layer{n}{1} = 'REGIONNAME2'; window_DV{n} = 1.95; window_AP{n} = 2.60;  n = n+1;
% add as many mice as you have


showROI_area = 'A1'; % show only A1 cells

% keep information
mouse_A1_F0s = mouse; date_A1_F0s = date; layer_A1_F0s = layer; window_DV_A1_F0s = window_DV; window_AP_A1_F0s = window_AP; n = 1; 



%% get roi_arealabel
ROIs_to_show_area = [];
for m = 1:length(date)
    for region = 1:length(layer{m})
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'labels.arearoilabel'],'-mat');
        roi_arealabels = roi_labels;
        
        % choose labels
        Allcells_area = (1:length(roi_arealabels))';
        nolabelcells_area = cellfun(@(x) isempty(x), roi_arealabels);
        A1cells = cellfun(@(x) any(strcmp('A1',x)), roi_arealabels);
        A2cells = cellfun(@(x) any(strcmp('A2',x)), roi_arealabels);
        
        switch showROI_area
            case 'All'
                ROIs_to_show_area{m}{region} = Allcells_area;
            case 'A1'
                ROIs_to_show_area{m}{region} = find(A1cells);
            case 'A2'
                ROIs_to_show_area{m}{region} = find(A2cells);
        end
    end
end


%% get roilabel
ROIs_to_show_celltype = [];
ROIs_to_show = [];
for m = 1:length(date)
    for region = 1:length(layer{m})
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'labels.roilabel'],'-mat');
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'baseline_F_level.mat']);
        
        Allcells = (1:length(roi_labels))';
        nolabelcells = cellfun(@(x) isempty(x), roi_labels);
        GADcells = cellfun(@(x) any(strcmp('GAD',x)),roi_labels);
        PVcells = cellfun(@(x) any(strcmp('PV',x)),roi_labels);
        SOMcells = cellfun(@(x) any(strcmp('SOM',x)),roi_labels);
        VIPcells = cellfun(@(x) any(strcmp('VIP',x)),roi_labels);
        
        filledcells = cellfun(@(x) any(strcmp('filled',x)),roi_labels);
        dimcells = cellfun(@(x) any(strcmp('reject_dim',x)),roi_labels);
        nancells = cellfun(@(x) any(strcmp('nan_ROI',x)),roi_labels);
        
        % set threshold for baseline F
        Fthreshold_dimcells = baselineF.F < F_threshold;
        
        rejectcells = filledcells|dimcells|nancells|Fthreshold_dimcells; % Can reject the GAD+ cells too
        
        
        switch showROI
            case 'All'
                ROIs_to_show_celltype{m}{region} = find(Allcells);
            case 'All_accept'
                ROIs_to_show_celltype{m}{region} = find(Allcells & ~rejectcells);
            case 'nolabel' % here, SCNN is included in nolabel cells
                ROIs_to_show_celltype{m}{region} = find(nolabelcells & ~rejectcells);
            case 'GAD'
                ROIs_to_show_celltype{m}{region} = find(GADcells & ~rejectcells);
            case 'PV'
                ROIs_to_show_celltype{m}{region} = find(PVcells & ~rejectcells);
            case 'SOM'
                ROIs_to_show_celltype{m}{region} = find(SOMcells & ~rejectcells);
            case 'VIP'
                ROIs_to_show_celltype{m}{region} = find(VIPcells & ~rejectcells);
        end
    end
end


%% overall ROIs to choose
for m = 1:length(date)
    for region = 1:length(layer{m})
        ROIs_to_show{m}{region} = intersect(ROIs_to_show_celltype{m}{region}, ROIs_to_show_area{m}{region});
    end
end



%% calculate area for each mouse, A1

area_all_A1 = [];

for m = 1:length(mouse)
    for region = 1:length(layer{m})

        [area_curr] = harmonicsF0_get_mean_responses(mouse{m}, date{m}, layer{m}{region});
        area_all_A1{m}{region} = area_curr;

    end
end


%% calculate coincidence-preference index and concatenate

coincidence_ind_A1_F0s = []; % coincidence_ind_A1{mouse}(ROIs)

for m = 1:length(area_all_A1)
    coincidence_ind_A1_F0s{m} = [];
    for region = 1:length(area_all_A1{m})
        
        % exc_detected file
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep 'harmonics_F0s\exc_detected_F.mat'],'response_ON');
        
        for ROI = 1:size(area_all_A1{m}{region},3)
            
            if ~ismember(ROI,ROIs_to_show{m}{region}) % rejected ROI. skip
                continue
            end
            
            response_ON_curr = response_ON{ROI}(:,3:end-1);
            % Input data are assumed to contain 7 delta-onset conditions
            % (-45:15:45 ms), plus two individual harmonic halves presented individually (first two columns).

            if sum(sum(response_ON_curr)) < area_threshold % not enough response number. skip
                continue
            end
            
            for F0_ind = 1:3 
                
                area_curr = area_all_A1{m}{region}(F0_ind,3:end-1,ROI); % exclude nosound trial
                
                
                % calculate CI individually and then average
                for dT_ind = 1:length(include_dTind)
                    dT_ind_curr = include_dTind(dT_ind);
                    if any(response_ON_curr(F0_ind,[4 dT_ind_curr])>0)
                        
                        if include_negative == 0
                            area_curr(area_curr<0) = 0;
                        end
                        if force_nonresponsive_to_zero == 1
                            area_curr = area_curr.*response_ON_curr(F0_ind,:);
                        end
                        
                        resp_coincident = max([0, area_curr(4)]);
                        resp_shifted = max([0, area_curr(dT_ind_curr)]);
                        
                        coincidence_ind_A1_F0s{m} = [coincidence_ind_A1_F0s{m};(resp_coincident-resp_shifted)/(resp_coincident+resp_shifted)];
                        
                    end
                end
            end
            
        end
    end
end




%% A2 data for harmonics_F0 experiments

mouse = []; date = []; layer = []; window_DV = []; window_AP = []; n = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%
% % A2 L2/3 pyramidal cells
% % extract sampling rate info from tiff header
% % %%%%%%%%%%%%%%%%%%%%%%%%%

mouse{n} = 'MOUSE3'; date{n} = 'DATE3'; layer{n}{1} = 'REGIONNAME3'; window_DV{n} = 1.61; window_AP{n} = 2.71;  n = n+1; 
mouse{n} = 'MOUSE4'; date{n} = 'DATE4'; layer{n}{1} = 'REGIONNAME4'; window_DV{n} = 1.95; window_AP{n} = 2.60;  n = n+1;
% add as many mice as you have



showROI_area = 'A2'; % show only the A2 cells

% keep information
mouse_A2_F0s = mouse; date_A2_F0s = date; layer_A2_F0s = layer; window_DV_A2_F0s = window_DV; window_AP_A2_F0s = window_AP; n= 1; 


%% get roi_arealabel
ROIs_to_show_area = [];
ROIs_to_show = [];
for m = 1:length(date)
    for region = 1:length(layer{m})
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'labels.arearoilabel'],'-mat');
        roi_arealabels = roi_labels;
        
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'baseline_F_level.mat']);
        
        % choose labels
        Allcells_area = (1:length(roi_arealabels))';
        nolabelcells_area = cellfun(@(x) isempty(x), roi_arealabels);
        A1cells = cellfun(@(x) any(strcmp('A1',x)), roi_arealabels);
        A2cells = cellfun(@(x) any(strcmp('A2',x)), roi_arealabels);
        
        switch showROI_area
            case 'All'
                ROIs_to_show_area{m}{region} = Allcells_area;
            case 'A1'
                ROIs_to_show_area{m}{region} = find(A1cells);
            case 'A2'
                ROIs_to_show_area{m}{region} = find(A2cells);
        end
    end
end

%% get label
ROIs_to_show_celltype = [];
for m = 1:length(date)
    for region = 1:length(layer{m})
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'labels.roilabel'],'-mat');
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep  'baseline_F_level.mat']);
        
        
        Allcells = (1:length(roi_labels))';
        nolabelcells = cellfun(@(x) isempty(x), roi_labels);
        GADcells = cellfun(@(x) any(strcmp('GAD',x)),roi_labels);
        PVcells = cellfun(@(x) any(strcmp('PV',x)),roi_labels);
        SOMcells = cellfun(@(x) any(strcmp('SOM',x)),roi_labels);
        VIPcells = cellfun(@(x) any(strcmp('VIP',x)),roi_labels);
        
        filledcells = cellfun(@(x) any(strcmp('filled',x)),roi_labels);
        dimcells = cellfun(@(x) any(strcmp('reject_dim',x)),roi_labels);
        nancells = cellfun(@(x) any(strcmp('nan_ROI',x)),roi_labels);
        
        
        % set threshold for baseline F
        Fthreshold_dimcells = baselineF.F < F_threshold;
        
        rejectcells = filledcells|dimcells|nancells|Fthreshold_dimcells;
        
        
        switch showROI
            case 'All'
                ROIs_to_show_celltype{m}{region} = find(Allcells);
            case 'All_accept'
                ROIs_to_show_celltype{m}{region} = find(Allcells & ~rejectcells);
            case 'nolabel' % here, SCNN is included in nolabel cells
                ROIs_to_show_celltype{m}{region} = find(nolabelcells & ~rejectcells);
            case 'GAD'
                ROIs_to_show_celltype{m}{region} = find(GADcells & ~rejectcells);
            case 'PV'
                ROIs_to_show_celltype{m}{region} = find(PVcells & ~rejectcells);
            case 'SOM'
                ROIs_to_show_celltype{m}{region} = find(SOMcells & ~rejectcells);
            case 'VIP'
                ROIs_to_show_celltype{m}{region} = find(VIPcells & ~rejectcells);
        end
    end
end


%% overall ROIs to choose
for m = 1:length(date)
    for region = 1:length(layer{m})
        ROIs_to_show{m}{region} = intersect(ROIs_to_show_celltype{m}{region}, ROIs_to_show_area{m}{region});
    end
end

%% calculate area for each mouse, A2

area_all_A2 = [];

for m = 1:length(mouse)
    for region = 1:length(layer{m})

        [area_curr] = harmonicsF0_get_mean_responses(mouse{m}, date{m}, layer{m}{region});
        area_all_A2{m}{region} = area_curr;
        
    end
end



%% calculate coincidence-preference index and concatenate

coincidence_ind_A2_F0s = []; % coincidence_ind_A1{mouse}(ROIs)

for m = 1:length(area_all_A2)
    coincidence_ind_A2_F0s{m} = []; 
    for region = 1:length(area_all_A2{m})
        
        % exc_detected file
        load([parentdir filesep date{m} filesep mouse{m} filesep layer{m}{region} filesep 'harmonics_F0s\exc_detected_F.mat'],'response_ON');
        
        for ROI = 1:size(area_all_A2{m}{region},3)
            
            if ~ismember(ROI,ROIs_to_show{m}{region}) % rejected ROI. skip
                continue
            end
            
            response_ON_curr = response_ON{ROI}(:,3:end-1);
            % Input data are assumed to contain 7 delta-onset conditions
            % (-45:15:45 ms), plus two individual harmonic halves presented individually (first two columns).
            
            if sum(sum(response_ON_curr)) < area_threshold % not enough response number. skip
                continue
            end
            
            for F0_ind = 1:3
                
                area_curr = area_all_A2{m}{region}(F0_ind,3:end-1,ROI); % exclude nosound trial
                
                % calculate CI individually and then average
                for dT_ind = 1:length(include_dTind)
                    dT_ind_curr = include_dTind(dT_ind);
                    if any(response_ON_curr(F0_ind,[4 dT_ind_curr])>0)
                        
                        if include_negative == 0
                            area_curr(area_curr<0) = 0;
                        end
                        if force_nonresponsive_to_zero == 1
                            area_curr = area_curr.*response_ON_curr(F0_ind,:);
                        end
                        
                        resp_coincident = max([0, area_curr(4)]);
                        resp_shifted = max([0, area_curr(dT_ind_curr)]);
                        
                        coincidence_ind_A2_F0s{m} = [coincidence_ind_A2_F0s{m};(resp_coincident-resp_shifted)/(resp_coincident+resp_shifted)];
                        
                    end
                end
            end
        end
    end
end


%% plot mouse vs window size

CI_A1_all = cellfun(@(x) nanmean(x), coincidence_ind_A1_F0s); % CI_A1(mouse)
CI_A2_all = cellfun(@(x) nanmean(x), coincidence_ind_A2_F0s); % CI_A2(mouse)

window_AP_A1_all = cell2mat(window_AP_A1_F0s);
window_DV_A1_all = cell2mat(window_DV_A1_F0s);
window_AP_A2_all = cell2mat(window_AP_A2_F0s);
window_DV_A2_all = cell2mat(window_DV_A2_F0s);

windowarea_A1_all = pi*(window_AP_A1_all/2).*(window_DV_A1_all/2);
windowarea_A2_all = pi*(window_AP_A2_all/2).*(window_DV_A2_all/2);

 
h = figure('position',[400 100 340 300]);
hold on
scatter(window_DV_A2_all, CI_A2_all,'markeredgecolor',[0 0.6 0.6]); 
xlabel('Window D-V size (mm)');
ylabel('Coincidence Pref Index');
pl_A2 = polyfit(window_DV_A2_all(~isnan(window_DV_A2_all)),CI_A2_all(~isnan(window_DV_A2_all)),1);
[R,p] = corrcoef(window_DV_A2_all(~isnan(window_DV_A2_all)),CI_A2_all(~isnan(window_DV_A2_all)));
line([1.4 2.8],[1.4 2.8]*pl_A2(1)+pl_A2(2),'color',[0 0.6 0.6]);
line([1.2 3],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([1.2 3]);
ylim([-0.4 0.4]);

h2 = figure('position',[700 100 400 300]);
hold on
scatter(window_AP_A2_all, CI_A2_all,'markeredgecolor',[0 0.6 0.6]); 
xlabel('Window A-P size (mm)');
ylabel('Coincidence Pref Index');
pl_A2 = polyfit(window_AP_A2_all(~isnan(window_AP_A2_all)),CI_A2_all(~isnan(window_AP_A2_all)),1);
[R,p] = corrcoef(window_AP_A2_all(~isnan(window_AP_A2_all)),CI_A2_all(~isnan(window_AP_A2_all)));
line([2.3 3.7],[2.3 3.7]*pl_A2(1)+pl_A2(2),'color',[0 0.6 0.6]);
line([2 4],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2 4]);ylim([-0.4 0.4]);

% use window area
h3 = figure('position',[700 100 340 300]);
hold on
scatter(windowarea_A2_all, CI_A2_all,'markeredgecolor',[0 0.6 0.6]); 
xlabel('Window area (mm2)');
ylabel('Coincidence Pref Index');
pl_A2 = polyfit(windowarea_A2_all(~isnan(windowarea_A2_all)),CI_A2_all(~isnan(windowarea_A2_all)),1);
[R,p] = corrcoef(windowarea_A2_all(~isnan(windowarea_A2_all)),CI_A2_all(~isnan(windowarea_A2_all)));
line([2.5 7.5],[2.5 7.5]*pl_A2(1)+pl_A2(2),'color',[0 0.6 0.6]);
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlim([2.5 7.5]);ylim([-0.4 0.4]);



%% alternatively, plot all cells

% DV size
h4 = figure('position',[100 400 900 300]);
subplot(1,2,1); hold on;
x = []; y = [];
for m = 1:length(mouse_A1_F0s)
    scatter(window_DV_A1_F0s{m}*ones(1,length(coincidence_ind_A1_F0s{m}))+0.025*randn(1,length(coincidence_ind_A1_F0s{m})),coincidence_ind_A1_F0s{m},'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor','none');
    x = [x, window_DV_A1_F0s{m}*ones(1,length(coincidence_ind_A1_F0s{m}))];
    y = [y, coincidence_ind_A1_F0s{m}'];
end
xlim([1.4 3]);
pl_A1_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
line([1.4 2.8],[1.4 2.8]*pl_A1_cells(1)+pl_A1_cells(2),'color','k');
line([1.4 2.8],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window D-V size (mm)');
ylabel('Coincidence Pref Index');
xlim([1 3]);

subplot(1,2,2); hold on;
x = []; y = [];
for m = 1:length(mouse_A2_F0s)
    scatter(window_DV_A2_F0s{m}*ones(1,length(coincidence_ind_A2_F0s{m}))+0.025*randn(1,length(coincidence_ind_A2_F0s{m})),coincidence_ind_A2_F0s{m},'markeredgecolor',[0 0.7 0.7],'markerfacecolor','none');
    x = [x, window_DV_A2_F0s{m}*ones(1,length(coincidence_ind_A2_F0s{m}))];
    y = [y, coincidence_ind_A2_F0s{m}'];
end
xlim([1.4 3]);
pl_A2_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
line([1.4 2.8],[1.4 2.8]*pl_A2_cells(1)+pl_A2_cells(2),'color','k');
line([1.4 2.8],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window D-V size (mm)');
ylabel('Coincidence Pref Index');
xlim([1 3]);

% AP size
h5 = figure('position',[100 700 900 300]);
subplot(1,2,1); hold on;
x = []; y = [];
for m = 1:length(mouse_A1_F0s)
    scatter(window_AP_A1_F0s{m}*ones(1,length(coincidence_ind_A1_F0s{m}))+0.025*randn(1,length(coincidence_ind_A1_F0s{m})),coincidence_ind_A1_F0s{m},'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor','none');
    x = [x, window_AP_A1_F0s{m}*ones(1,length(coincidence_ind_A1_F0s{m}))];
    y = [y, coincidence_ind_A1_F0s{m}'];
end
xlim([2 4]);
pl_A1_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
line([2.3 3.7],[2.3 3.7]*pl_A1_cells(1)+pl_A1_cells(2),'color','k');
line([2 4],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window A-P size (mm)');
ylabel('Coincidence Pref Index');

subplot(1,2,2); hold on;
x = []; y = [];
for m = 1:length(mouse_A2_F0s)
    scatter(window_AP_A2_F0s{m}*ones(1,length(coincidence_ind_A2_F0s{m}))+0.025*randn(1,length(coincidence_ind_A2_F0s{m})),coincidence_ind_A2_F0s{m},'markeredgecolor',[0 0.7 0.7],'markerfacecolor','none');
    x = [x, window_AP_A2_F0s{m}*ones(1,length(coincidence_ind_A2_F0s{m}))];
    y = [y, coincidence_ind_A2_F0s{m}'];
end
xlim([2 4]);
pl_A2_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
line([2.3 3.7],[2.3 3.7]*pl_A2_cells(1)+pl_A2_cells(2),'color','k');
line([2 4],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window A-P size (mm)');
ylabel('Coincidence Pref Index');



% plot against window area

windowarea_A1_F0s = [];
for m = 1:length(window_AP_A1_F0s)
    windowarea_A1_F0s(m)= pi*(window_DV_A1_F0s{m}/2)*(window_AP_A1_F0s{m}/2);
end
windowarea_A2_F0s = [];
for m = 1:length(window_AP_A2_F0s)
    windowarea_A2_F0s(m)= pi*(window_DV_A2_F0s{m}/2)*(window_AP_A2_F0s{m}/2);
end

h6 = figure('position',[100 1000 850 300]);
subplot(1,2,1); hold on;
x = []; y = [];
for m = 1:length(mouse_A1_F0s)
    scatter(windowarea_A1_F0s(m)*ones(1,length(coincidence_ind_A1_F0s{m}))+0.1*randn(1,length(coincidence_ind_A1_F0s{m})),coincidence_ind_A1_F0s{m},'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor','none');
    x = [x, windowarea_A1_F0s(m)*ones(1,length(coincidence_ind_A1_F0s{m}))];
    y = [y, coincidence_ind_A1_F0s{m}'];
end
xlim([2.5 7.5]);
pl_A1_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
[R p] = corrcoef(x(~isnan(x+y)),y(~isnan(x+y)));
line([2.5 7.5],[2.5 7.5]*pl_A1_cells(1)+pl_A1_cells(2),'color','k');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window area (mm2)');
ylabel('Coincidence Pref Index');

subplot(1,2,2); hold on;
x = []; y = [];
for m = 1:length(mouse_A2_F0s)
    scatter(windowarea_A2_F0s(m)*ones(1,length(coincidence_ind_A2_F0s{m}))+0.1*randn(1,length(coincidence_ind_A2_F0s{m})),coincidence_ind_A2_F0s{m},'markeredgecolor',[0 0.7 0.7],'markerfacecolor','none');
    x = [x, windowarea_A2_F0s(m)*ones(1,length(coincidence_ind_A2_F0s{m}))];
    y = [y, coincidence_ind_A2_F0s{m}'];
end
xlim([2.5 7.5]);
pl_A2_cells = polyfit(x(~isnan(x+y)),y(~isnan(x+y)),1);
[R p] = corrcoef(x(~isnan(x+y)),y(~isnan(x+y)));
line([2.5 7.5],[2.5 7.5]*pl_A2_cells(1)+pl_A2_cells(2),'color','k');
line([2.5 7.5],[0 0],'linestyle',':','color',[0.5 0.5 0.5]);
xlabel('Window area (mm2)');
ylabel('Coincidence Pref Index');
