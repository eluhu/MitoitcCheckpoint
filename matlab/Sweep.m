function [ Results ] = Sweep( modelName, p_laserMax, p_laserMin, k_tensRange, p_cutXRange )
%SWEEP Is a fucntion that explores the parameter space for laser cutting
%experiments in a systematic way. It takes the model name ('modelName)',
%the maximum time at which the laser cuts ('p_laser_ini'),
%a vector holding the set of Xsome reattachment rates (k_tensRange) to be considered
%and a vector holding the set of degrees of Xsome detachment to be considered ('p_cutXRange') as arguments
%Returns the timecourses as well as summary datasets for the timepoints at
%which the checkpoint becomes no longer responsive to reengagement
%   Detailed explanation goes here

if k_tensRange == 0
    k_tensRange = 0.00:0.005:0.1;
end

if p_cutXRange == 0
    p_cutXRange = 0:0.05:0.975;
end

%Setup routine
n = 1;
idx = 1;
p_laser = p_laserMax;
tog = 1;                        %toggle indicating first finding: tog=1 no parameter set in which CP is responsive to re-rengagement has been found yet, tog=0 set found
model = modelName;              %name of the model: HEslow, HEfast, HEnpp, SIMM, SIMMpp, SIMMmm
lookingat = 'ResToReengagment'; %i.e. Responsiveness to Reengagement
today = datestr(clock,29);      %Converts todays date into string of form yyyy-mm-dd
size = length(k_tensRange) * length(p_cutXRange);
modelHandle = str2func(model);  %Converts string identifying model into function handle used to call model function
Results = zeros(size, 4);       %preallocates memory for the results


%Sweep over rate of Xsome reattachment
for k_tens = k_tensRange
    
    %Sweep over extent of Xsomes cut by laser (cutX = fraction of Xsomes attached after cutting)
    for p_cutX = p_cutXRange
        
        %As long as no conditions under which CP disengagement is
        %reversible have been found (tog = 1)...
        while tog == 1
             %...solve the model corresponging to modelName for the given
             %set of parameters and re-evaluate whether CP is responsive
             [tY, op, d1MCCt] = modelHandle(p_laser, p_cutX, k_tens, n);
             
             %if the model returns op=1 (checkpoint is responsive!) set
             % write timecoure to [...]-data.dat file and Results-entry to [...]_summary.dat file
             % and exit while loop
             if op == 1 && tog == 1
                 Results(idx,:) = [n, k_tens, p_cutX, p_laser];
                 tog = 0;
                 filenameData = sprintf('%s_%s_%s_#%d_data.dat',today, model, lookingat, n);
                 save(filenameData, 'tY', '-ASCII');
                 filenameSummary = sprintf('%s_%s_%s_#%d_summary.dat',today, model, lookingat, n);
                 save(filenameSummary, 'Results', '-ASCII');
             %if the tested time of laser cutting is below a set minimum,
             %but the checkpoint is not responsive, do the same as above.
             elseif op == 0 && tog == 1 && p_laser < p_laserMin
                 Results(idx,:) = [n, k_tens, p_cutX, 20];
                 tog = 0;
                 filenameData = sprintf('%s_%s_%s_#%d_data.dat',today, model, lookingat, n);
                 save(filenameData, 'tY', '-ASCII');
                 filenameSummary = sprintf('%s_%s_%s_#%d_summary.dat',today, model, lookingat, n);
                 save(filenameSummary, 'Results', '-ASCII');
             end
             %otherwise stay in the loop and cut one time increment
             %earlier
             n = n + 1;
             p_laser = p_laser - 0.1;
        end
        
        %reset the mechanism and set counter for next p_cutX
        idx = idx + 1;
        tog = 1;
    end
    %...after all p_cutX have been tested for a single k_tens, reset
    %p_laser to its maximum and go to next k_tens
    p_laser = p_laserMax;
end

%--------------------------------------------------------------------------------------------------------
% SUMMARY
%--------------------------------------------------------------------------------------------------------
i = 1;
summary = zeros(size, 10);

%for each pair of p_cutX and k_tens... 
for i = 1:size
    
    %...open the corresponding timecourse file
    filename = sprintf('%s_%s_%s_#%d_data.dat',today, model, lookingat, Results(i,1));
    C = importdata(filename);
    %...find the index at which Xtens has been reduced to less than 1
    notOne = find(C(:,8)~=1);
    changeID = min(notOne);
    %...and write the conditions of the system at that timepoint into the
    %summary array
    summary(i,1) = Results(i,2);
    summary(i,2:10) = C(changeID,:);

end

%save summaryarray
summaryname = sprintf('%s_ReengRespo_SummaryData.txt', model);
save(summaryname, 'summary', '-ASCII');

end

