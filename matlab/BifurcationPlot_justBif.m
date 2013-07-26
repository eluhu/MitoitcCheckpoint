function [ ] = BifurcationPlot_justBif( model )
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
    conditions = {'CycB047(SS)', 'CycB030', 'CycB020', 'CycB010', 'CycB008', 'CycB006', 'CycB004', 'CycBvar'};    
elseif modelID == 5
    conditions = {'CycB077(SS)', 'CycB040', 'CycB020', 'CycB015', 'CycB011', 'CycB005', 'CycB002', 'CycBvar'};
elseif modelID == 6 || modelID == 7
    conditions = {'CycB047(SS)', 'CycB020', 'CycB010', 'CycB008', 'CycB006', 'CycB004', 'CycB0028' 'CycBvar'};
else
    conditions = {'CycB077(SS)', 'CycB040', 'CycB030', 'CycB020', 'CycB015', 'CycB011', 'CycB010', 'CycBvar'};
end

% Define a set of 7 shades of grey and one dark red, which will be used to
% colour the bifurcation diagrams displayed in the background

greys = {[0.4 0.4 0.4], [0.5 0.5 0.5], [0.6 0.6 0.6], [0.7 0.7 0.7], [0.8 0.8 0.8],[0.9 0.9 0.9], [0.95 0.95 0.95], [0.42 0 0]};

modelColor = {[0 0.8039 0.4000], [0.1098 0.5254 0.9333], [0.6980 0.1333 0.1333], [1.0000 0.3804 0.0118], [0.5764 0.4392 0.8588], [0.42 0 0], [0.8039 0.1608 0.5647]};
greys{8} = modelColor{modelID};


%Define the dimensions of the figure
set(figure(1),'Position', [100 100 1000 600]);
figure(1);

%Loop over conditions
for i = 1:8
    
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
    
    hold on
    
    %Plot the upper stable steady state...
    h(i,1) = plot(C(1:unIDmin-1,1),C(1:unIDmin-1,2));
    %... the unstable steady state (as dashed line)...
    h(i,2) = plot(C(unIDmin:unIDmax,1),C(unIDmin:unIDmax,2),'--');
    %... and the lower stable steady state.
    h(i,3) = plot(C(unIDmax+1:length(C),1),C(unIDmax+1:length(C),2));
    
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

%LABELS, AXIS RANGES AND TITLE
xlabel(parX);
ylabel(varY);
xlim([0 1]);
ylim([0 2]);
plottitle = sprintf('%s',model);
title(plottitle,...
    'FontWeight','bold')

set(findall(figure(1),'type','text'),'FontName', 'Rockwell');

end



