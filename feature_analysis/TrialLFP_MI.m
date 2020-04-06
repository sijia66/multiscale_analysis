%use MI to calculate the distribution
%
%first, we need to bin this into boxes. They only accept discrete inputs

X = ceil(s_z_targ_ML_window);
Y = s_z_targ;
mi_x_y = mi(X,Y);    
drive_idx = 1;

depthProfile = (Trials(1).Depth{1,drive_idx})'; % in micron
loweredElectrodes  = find(depthProfile > 0);

mi_x_y = zeros(size(X,2),size(X,3));

for ei = 1:size(X,2)
    for fi = 1:size(X,3)
        
        if (drive_idx == 1)
            
            datTemp = squeeze(X(:,ei,fi));
            mi_x_y(ei,fi) = mi(datTemp,Y);
            
        elseif (drive_idx == 2)
            if(ismember(ei,loweredElectrodes))
                datTemp = squeeze(X(:,ei,fi));
                mi_x_y(ei,fi) = mi(datTemp,Y);
            else
                mi_x_y(ei,fi) = 0;
            end
        else
            disp('Check your drive index')
            break
        end
    end
end
disp('All processing done!')

figure
imagesc(mi_x_y')
colorbar
axis xy
xlabel('Electrode Number')
ylabel('Frequency (Hz)')

yticklabels = 0:5:300;
yticks = linspace(1, size(mi_x_y, 2), numel(yticklabels));
set(gca, 'yTick', yticks, 'yTickLabel', yticklabels)


%%
%choose a frequency band and plot as a bubble plot

%%
% we are going to calculate the mutual information between two drives
% first, we can do this for a particular frequency band and 

perIndex = 95;

mi_x_y = calMI_xy(s_z_targ_ECOG',s_z_targ_SC32',trialInfo.badECOG, trialInfo.badSC32);
imagesc(mi_x_y')
xlabel('ECOG electrodes')
ylabel('SC32 electrodes')
colorbar


cohSeq_fr_temp = mi_x_y;
cohSeq_fr_temp(31:32,:) = 0;
[nRow,nCol] = size(cohSeq_fr_temp);
cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol,1);

coh_th(3) = prctile(cohSeq_fr_temp, perIndex );
fprintf('ECOG vs. SC32: %.4f\n',coh_th(3));

cMatrix = mi_x_y';
cMatrix(31:32,:) = 0; % this is an ad hoc solutions
graphParam.threshold = coh_th(3);
graphParam.thresholdPercentile = perIndex;

graphParam.driveName1 = 'SC32';
graphParam.P_chamber1 = P_chamber_SC32;

graphParam.driveName2 = 'ECOG';
graphParam.P_chamber2 = P_chamber_ECOG_3;

graphParam.freqRange = freqRange;
graphParam.labelElectrode = true;
graphParam.drive_idx =  2; %1 for ECOG
graphParam.textShift = 0.1;


figure
A = cMatrix > graphParam.threshold;
inter_gplot(A,graphParam)








