function  mapCoherenceGraph(cMatrix,P_chamber, graphParam)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

th = graphParam.threshold;
textShift =  graphParam.textShift;
driveName = graphParam.driveName;
drive_idx = graphParam.drive_idx;
labelElectrode = graphParam.labelElectrode;
acList = graphParam.acList;
freqRange = graphParam.freqRange;

%let's apply a threshold
cMatrix(isnan(cMatrix)) = 0;
cMatrix(cMatrix > th) = 1;
cMatrix(cMatrix ~= 1) = 0;

%we need to complete the coherence matrix
cMatrix = cMatrix + cMatrix';

% plot and annotate
gplot(cMatrix,P_chamber,'-*')
axis square 
axis([-6 6 -6 6])
xlabel('mm')
ylabel('mm')
txt = sprintf('Cutoff Threshold %0.2f (%d th Percentile)',graphParam.threshold,graphParam.thresholdPercentile);
text(-3,-8,txt,'FontSize',14)


txt = sprintf('Averaged Between %0.2f and %0.2f Hz',freqRange(1),freqRange(2));
text(-3,-8.8,txt,'FontSize',14)


%label the electrode
if labelElectrode
    hold on
   
    %need to fix this for SC32
    for i = 1:length(acList)
        text(P_chamber(i,1)+textShift,P_chamber(i,2),num2str(i))
    end
    hold off

end

end

