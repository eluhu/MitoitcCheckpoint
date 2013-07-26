function [ ] = BifurcationPlot3D( model )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Clear open figures
figure(1);
clf;

%Define what's going on... stings will be used to handle file input
%(filenames) and plot titles as well as axis lables
what = 'Bifurcation';
parX = 'Xtens';
varY = 'MCCt';

models = {'HEslow', 'HEfast', 'SIMM', 'SIMMmm', 'HEnpp', 'SIMMRegPP', 'SIMMpp'};

%Compare the input argument 'model' to the list of supported models
%specified in "models" to pick the appropriate colormap
for k = 1:7
    if strcmp(model, models{k}) == 1
        modelID = k;
    end
end

%Specify the colormaps:
%He model: a set of six greens (dark to light)
colormap{1} = {[0 0.39 0], [0.18 0.54 0.34], [0.23 0.70 0.44], [0.40 0.80 0.66], [0.6078 0.8039 0.6078], [0.60 0.98 0.60] };
%HeFast model: a set of six blues (dark to light)
colormap{2} = {[0.1 0.1 0.44], [0.27 0.5 0.7], [0.25 0.41 0.88], [0.39 0.58 0.93], [0.68 0.76 0.87], [0.68 0.84 0.90]};
%SIMM model: a set of six reds (dark to light)
colormap{3} = {[0.52 0 0], [0.70 0.13 0.13], [0.86 0.07 0.23], [0.94 0.5 0.5] [0.98 0.5 0.44] [1 0.63 0.48]};
%SIMMhelena model: a set of six yellows/browns (dark to light)
colormap{4} = {[0.82 0.41 0.12], [1 0.54 0], [1 0.64 0], [1 0.84 0], [1 0.90 0.71], [1 1 0.87] };
%He model: a set of six violets (dark to light)
colormap{5} = {[0.2941 0 0.5098], [0.4078 0.0941 0.5451], [0.6039 0.1961 0.8039], [0.6980 0.2275 0.9333], [0.7490 0.2431 1.0000], [0.8784 0.4000 1.0000]};
%SIMMRegPP model: a set of six violets (dark to light)
colormap{6} = {[0.2941 0 0.5098], [0.4078 0.0941 0.5451], [0.6039 0.1961 0.8039], [0.6980 0.2275 0.9333], [0.7490 0.2431 1.0000], [0.8784 0.4000 1.0000]}'';
%SIMMRegPP model: a set of six violets (dark to light)
colormap{7} = {[0.2941 0 0.5098], [0.4078 0.0941 0.5451], [0.6039 0.1961 0.8039], [0.6980 0.2275 0.9333], [0.7490 0.2431 1.0000], [0.8784 0.4000 1.0000]}'';

%Find the appropriate set of conditions used to generate the bifurcation
%diagrams (model dependent) - stcmp basically just compares the strings
if modelID == 3 || modelID == 4
    conditions = {'CycB047(SS)', 'CycB030', 'CycB020', 'CycB010', 'CycB008', 'CycB006', 'CycB004'};
    cyclinLevels = [0.47, 0.30, 0.20, 0.10, 0.08, 0.06, 0.04];
elseif modelID == 5
    conditions = {'CycB077(SS)', 'CycB040', 'CycB020', 'CycB015', 'CycB011', 'CycB005', 'CycB002'};
    cyclinLevels = [0.77, 0.40, 0.20, 0.15, 0.11, 0.05, 0.002];
elseif modelID == 6 || modelID == 7
    conditions = {'CycB047(SS)', 'CycB020', 'CycB010', 'CycB008', 'CycB006', 'CycB004', 'CycB0028'};
    cyclinLevels = [0.47, 0.20, 0.10, 0.08, 0.06, 0.04, 0.028];
else
    conditions = {'CycB077(SS)', 'CycB040', 'CycB030', 'CycB020', 'CycB015', 'CycB011', 'CycB010'};
    cyclinLevels = [0.77, 0.40, 0.30, 0.20, 0.15, 0.11, 0.10];
end

% Define a set of 7 shades of grey and one dark red, which will be used to
% colour the bifurcation diagrams displayed in the background
greys = {[0.4 0.4 0.4], [0.5 0.5 0.5], [0.6 0.6 0.6], [0.7 0.7 0.7], [0.8 0.8 0.8],[0.9 0.9 0.9], [0.95 0.95 0.95], [0.42 0 0]};

%Define the dimensions of the figure
set(figure(1),'Position', [100 100 1000 600]);
figure(1);
 



%Loop over conditions
for i = 1:7
    
    %Build the input filename using the strings defined above and import
    %the file
    filename = sprintf('%s_%s_%s-%s_%s.dat', what, model, parX, varY, conditions{i});
    C = importdata(filename);
    %set the state identifier for the first entry to one 9otherwise it is always 2 =
    %unstable)
    C(1,4) = 1;
    %Using the state identifier (1 = stable, 2 = unstable), find the
    %indices of the first and last entry for the unstable branch
    unID = find(C(:,4)~=1);
    unIDmin = min(unID);
    unIDmax = max(unID);
    
    cyclinVector = ones(length(C),1) * cyclinLevels(i);
    
    
    hold on
    
    %Plot the upper stable steady state...
    h(i,1) = plot3(C(1:unIDmin-1,1),C(1:unIDmin-1,2), cyclinVector(1:unIDmin-1,1));
    %... the unstable steady state (as dashed line)...
    h(i,2) = plot3(C(unIDmin:unIDmax,1),C(unIDmin:unIDmax,2), cyclinVector(unIDmin:unIDmax,1));
    set(h(i,2),'LineStyle',':')
    %... and the lower stable steady state.
    h(i,3) = plot3(C(unIDmax+1:length(C),1),C(unIDmax+1:length(C),2),cyclinVector(unIDmax+1:length(C),1));
    
    
    
    %SETTING UP THE REPRESENTATION IN THE FIGURE LEGEND
    %For each condition defined in cell arrray 'conditions' group the plots
    %of the three branches and set this group to be the "parent" of the three subbranches 
    %and set the parent to be represented in the legend rather than each of
    %the branches individually
    bifurcation(i) = hggroup;
    set(h(i,1:3), 'Parent', bifurcation(i))
    set(get(get(bifurcation(i), 'Annotation'), 'LegendInformation'),'IconDisplayStyle','on');
    %Assing the colour according to the set of colour defined in 'greys'
    set(h(i,:), 'Color', greys{i}, 'LineWidth', 2);
    
    hold off
end

hold on

filename = sprintf('%s_%s_%s-%s_%s.dat', what, model,'Xtens', 'MCCt', 'CycBvar');
mcc_xtens = importdata(filename);
filename = sprintf('%s_%s_%s-%s_%s.dat', what, model, 'Xtens', 'CycB', 'MCCtvar');
cycB_xtens = importdata(filename);

mcc_xtens(1,4) = 1;
cycB_xtens(1,4) = 1;
%Using the state identifier (1 = stable, 2 = unstable), find the
%indices of the first and last entry for the unstable branch
unID_mcc = find(mcc_xtens(:,4)~=1);
unIDmin_mcc = min(unID_mcc);
unIDmax_mcc = max(unID_mcc);

unID_cycB = find(cycB_xtens(:,4)~=1);
unIDmin_cycB = min(unID_cycB);
unIDmax_cycB = max(unID_cycB);

if mcc_xtens(:,1) == cycB_xtens(:,1)

    xtens = mcc_xtens(:,1);
    mcc= mcc_xtens(:,2);
    cycB = cycB_xtens(:,2);
    
else
    fprintf('Data for 3D-bifurcation plot does not match')
end  
    
var(1) = plot3(xtens(1:unIDmin_mcc-1),mcc(1:unIDmin_mcc-1), cycB(1:unIDmin_cycB-1));
%... the unstable steady state (as dashed line)...
var(2) = plot3(xtens(unIDmin_mcc:unIDmax_mcc),mcc(unIDmin_mcc:unIDmax_mcc), cycB(unIDmin_cycB:unIDmax_cycB));
set(var(2),'LineStyle',':')
%... and the lower stable steady state.
var(3) = plot3(xtens(unIDmax_mcc+1:length(xtens)), mcc(unIDmax_mcc+1:length(mcc)), cycB(unIDmax_cycB+1:length(cycB)));

 variationBif = hggroup;
    set(var(1:3), 'Parent', variationBif)
    set(get(get(variationBif, 'Annotation'), 'LegendInformation'),'IconDisplayStyle','on');
    %Assing the colour according to the set of colour defined in 'greys'
    set(var(1:3), 'Color', greys{8}, 'LineWidth', 2);
%--------------------------------------------------------------------------
% Add Simulation Data
%--------------------------------------------------------------------------
hold on

%Build the input filename using the strings defined above and import
%the file
filename = sprintf('%s_ReengRespo_SummaryData.txt', model);
C = importdata(filename);
%Specify the set of k_tens conditions to be considered in the plot
plot_k_tens = [0, 0.02, 0.04, 0.06, 0.08, 0.1];

%Loop over each element to extract the value of each entry and generate the
%legend label for each series
for u = 1:length(plot_k_tens)
    legendDataP{u} = sprintf('k_{tens}: %.2f', plot_k_tens(u));
end

%Preallocate memory for the arraz holding the plotting data
datapoints = zeros(20,length(plot_k_tens));
maxSensPoint = zeros(1, length(plot_k_tens));
maxSensXtens = zeros(1, length(plot_k_tens));
maxSensMCCt = zeros(1, length(plot_k_tens));

%Loop over the conditions specified in p_k_tens
for m = 1:length(plot_k_tens)
    
    %Find the indices correspondinng to datapoints for the current k_tens
    %value
    plotrange = find(C(:,1) == plot_k_tens(m));
    
    %If the checkpoint disengagement is irreversible under the given
    %conditions, exclude points. (Judged based on threshold time set in Sweep: p_lasermin 3min for HeFast, 25min else)
    effectiveLength = length(plotrange);
    
    %If not HeFast
    if modelID ~= 2
        for p = 1:length(plotrange)
            %...if time determined for irrev < 25.5
            if C(plotrange(p), 2) < 25.5
                effectiveLength = effectiveLength - 1;
            end
        end
    %If it is HeFast
    else
        for p = 1:length(plotrange)
            %...if time determined for irrev <3
            if C(plotrange(p), 2) < 3
                effectiveLength = effectiveLength - 1;
            end
        end
    end
    
    %... reduce range to be plotted according to the determined effective
    %length
    plotrange = plotrange(1:effectiveLength);
        
    
    %Plot the MCCt (column 5) vs Xtens (column 9) in the specified range
    %marked by circles (no line)
    datapoints(:,m) = scatter3(C(plotrange,9), C(plotrange,5),C(plotrange,3),'o', 'MarkerFaceColor', colormap{modelID}{m}, 'MarkerEdgeColor', [0.3 0.3 0.3]);
    
    %Find the maximum sensitivity for reengagement
    maxSensXtens(m) = max(C(plotrange,9));
    maxSensMCCt(m) = max(C(plotrange,5));
    maxSensCycB(m) = max(C(plotrange,3));
    maxSensPoint(m) = scatter3(maxSensXtens(m), maxSensMCCt(m), maxSensCycB(m), 's', 'MarkerFaceColor', colormap{modelID}{m}, 'MarkerEdgeColor', [0.3 0.3 0.3]);
    
    dataSeries(m) = hggroup;
    set(datapoints(:,m), 'Parent', dataSeries(m));
    set(maxSensPoint(m), 'Parent', dataSeries(m));
    set(get(get(dataSeries(m), 'Annotation'), 'LegendInformation'),'IconDisplayStyle','on');
    
end


%Concatenate the cell arrays holding the names for the legend horizontally
legendNames = horzcat(conditions, legendDataP);

%Setup legend (names to be displayed, appearance, position)
legendHandle1 = legend(legendNames);
set(legendHandle1, 'Box', 'off');
set(legendHandle1, 'FontSize', 10);
set(legendHandle1, 'Location', 'EastOutside');

hold off

%LABELS, AXIS RANGES AND TITLE
xlabel(parX);
ylabel(varY);
zlabel('CycB');
set(gca,'ZDir','reverse', 'Zlim', [0 (max(cyclinLevels)+0.2)], 'Xlim', [0 1.05], 'Ylim', [0 2.05])
plottitle = sprintf('%s-Model 3D-%s Diagram: %s vs. %s vs. CycB', model, what, parX, varY);
title(plottitle,...
    'FontWeight','bold')

set(findall(figure(1),'type','text'),'FontName', 'Rockwell');

end



