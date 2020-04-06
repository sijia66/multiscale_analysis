data = dataTarget_SC32.target1;
dataTrials = dataTarget_SC32.target1_trials;
timeReachStart = 2000;
timeInterest = 300;

%normalize using timing between go to TargOn
ti = 1;
ei = 20;

freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParam.dn = 0.02;
freqParams.contFlag = 1;

freqRange = [0 50];   

dataTemp = squeeze(data(ti,ei,:));
%first normalize
dataTemp = dataTemp / norm(dataTemp);

timeAcqOn = dataTrials(ti).StartAq - dataTrials(ti).ReachStart + timeReachStart;
timeTargsOn = dataTrials(ti).TargsOn - dataTrials(ti).ReachStart + timeReachStart;
dataTemp = dataTemp(timeAcqOn:timeTargsOn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do frequency calculation

dat = dataTemp';
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
freqRangeIndex = cal_index_freq(freq ,freqRange(1),freqRange(2));

s_n_mean = squeeze(mean(s));
s_n_std = squeeze(std(s));

figure
plot(freq(freqRangeIndex ),10*log(s(:,freqRangeIndex)') )
xlabel('Frequency')
ylabel('Power Amplitude (dB)')
title('Frequency Spectra Between Go and TargsOn')

figure
imagesc(10*log(s(:,freqRangeIndex)'))
colorbar
axis xy

%%
% we are gonna get 
dataTemp = squeeze(data(ti,ei,:));
%first normalize
dataTemp = dataTemp / norm(dataTemp);


figure
plot(dataTemp)
dat = dataTemp';
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);

%get z-score
sN = size(s,1);
s_z = (s - repmat(s_n_mean,sN,1)) ./ repmat(s_n_std,sN,1);


figure
imagesc(s_z')
caxis([-10 10])
colorbar
axis xy


%%
% need to process the trials

timeReachStart = 2000;

freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParam.dn = 0.02;        
freqParams.contFlag = 1;
freqRange = [0 300];

% extract features
preMovementTime = 200; % in ms;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s_z_targ_features = [];
s_z_targ_ML = [];
s_z_targ = [];

for targi = 1:7
    
    eval(sprintf('data = dataTarget.target%d;',targi))
    eval(sprintf('dataTrials = dataTarget.target%d_trials;',targi))
    
    s_z_pool = [];
    s_z_targ_trial_features =  [];
    for ei = 1:211
        for ti = 1:size(data,1)
            
            dataTemp = squeeze(data(ti,ei,:));
            %first normalize
            dataTemp = dataTemp / norm(dataTemp);
            
            timeAqGo = dataTrials(ti).StartAq - dataTrials(ti).ReachStart + timeReachStart;
            timeGo = dataTrials(ti).Go - dataTrials(ti).ReachStart + timeReachStart; %relative to reach start
            timeTargsOn = dataTrials(ti).TargsOn - dataTrials(ti).ReachStart + timeReachStart;
            
            dataTemp = dataTemp(timeAqGo:timeTargsOn);
            
            %do frequency calculation
            dat = dataTemp';
            [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
            freqRangeIndex = cal_index_freq(freq ,freqRange(1),freqRange(2));
            
            s_n_mean = squeeze(mean(s));
            s_n_std = squeeze(std(s));
            
            % now look at the trials
            % we are gonna get
            dataTemp = squeeze(data(ti,ei,:));
            %(timeReachStart - timeInterest  + 1) : (timeReachStart + timeInterest) ));
            %first normalize
            dataTemp = dataTemp / norm(dataTemp);
            dat = dataTemp';
            [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
            
            %get z-score
            sN = size(s,1);
            s_z = (s - repmat(s_n_mean,sN,1)) ./ repmat(s_n_std,sN,1);
            
            %extract features from these plots
            s_length = size(s_z,1);
            s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.tapers(1) * 0.1;
            time_indices = find((s_timePoints <= timeReachStart / 1000) ...
                & (s_timePoints >= (timeReachStart - preMovementTime) / 1000)) ;
            freq_indices = find(freq <= freqRange(2)) ;
            s_z_preMovement = squeeze(mean(s_z(time_indices,freq_indices)));
            
            %s_z_targ_data(targi,:,:) = s_z_pool_mean;
            s_z_targ_trial_features(ti,ei,:) = s_z_preMovement;
  
        end
    end
    
    %build ML datasets
    s_z_targ_ML = cat(1,s_z_targ_ML,s_z_targ_trial_features);
    s_z_targ = [s_z_targ;ones(size(data,1),1) * targi];
    
    
    s_z_targ_features(targi,:,:) = squeeze(mean(s_z_targ_trial_features));
   
    disp(strcat('finished:',num2str(targi)))
    
end
disp('All Processing Done')

%%
%format into ML format 



%%
% do some plots

freqRange_display_pool = [0.5 3; 3 8; 8 12;12 38;38 42];
%freqRange_display_pool = [70 85];
labelElectrode = 0;
onlyLowered = 0;%0 for ECOG

acList = [experiment.hardware.microdrive(drive_idx).electrodes.acquisitionid];
magFac = 100;
textShift = 0.1;

%load the correct P_chamber index
if drive_idx == 1
    P_chamber = P_chamber_ECOG_3;
elseif drive_idx == 2
     P_chamber = P_chamber_SC32;
else
    disp('Check your index')
end


depthProfile = (Trials(1).Depth{1,drive_idx})'; % in micron
loweredElectrodes  = find(depthProfile > 0);
%then we do some plotting
     
for fri = 1:size(freqRange_display_pool,1)
    
    freqRange_display = freqRange_display_pool(fri,:);
    freqRange = find(freq <= freqRange_display(2) & freq >= freqRange_display(1));

    
    figure('units','normalized','outerposition',[0 0 1 1])
    for targi = 1:7
        

        
        subplot(2,4,targi)
        bubbleValues = squeeze(mean(s_z_targ_features(targi,:,freqRange),3));
        % process the bubble values, red as increase and blue means decrease
        %increasing
        
        if onlyLowered
            bubbleValues_inc_indices = intersect(find(bubbleValues >0), loweredElectrodes);
            bubbleValues_dec_indices = intersect(find(bubbleValues  < 0), loweredElectrodes);
        else
            bubbleValues_inc_indices = find(bubbleValues >0 );
            bubbleValues_dec_indices = find(bubbleValues < 0);
        end
        
        bubbleValues_zero_indices = find(bubbleValues == 0);
        
        scatter(P_chamber(bubbleValues_inc_indices,1),P_chamber(bubbleValues_inc_indices,2)...
            ,bubbleValues(bubbleValues_inc_indices) * magFac, ...
            'MarkerFaceColor','r')
        axis equal
        
        %label the negative z-scores
        hold on
        
        scatter(P_chamber(bubbleValues_dec_indices,1),P_chamber(bubbleValues_dec_indices,2)...
            ,abs(bubbleValues(bubbleValues_dec_indices)) * magFac,...
            'MarkerFaceColor','b','MarkerEdgeColor','none')
        

        scatter(P_chamber(bubbleValues_zero_indices,1),P_chamber(bubbleValues_zero_indices,2)...
            ,magFac  / 2, 'MarkerFaceColor','g')
        hold off
        title(sprintf('Target %d: Frequency Range %0.1f and %0.1f', ...
            targi,freqRange_display(1),freqRange_display(2)))
        xlabel('mm')
        ylabel('mm')
        
        
        %label the electrodes
        if labelElectrode
            
            hold on
            if drive_idx == 1 %ECOG 
                for i = 1:size(s_z_targ_features,2)
                    text(P_chamber(i,1)+textShift,...
                        P_chamber(i,2)....
                        ,num2str(i))
                end
            elseif drive_idx == 2 %SC32
                for i = 1:length(loweredElectrodes)
                    text(P_chamber(loweredElectrodes(i),1)+textShift,...
                        P_chamber(loweredElectrodes(i),2)....
                        ,num2str(i))
                end
            else
                disp('Check drive index number')
            end
            hold off
        end
    end
end
%%
figure
plotTF_dB(s_z_pool_mean(:,freqRangeIndex),freq,freqRange,freqParams.tapers(1))
caxis([-3 3])
colorbar
axis xy

title(sprintf('Target %d for Electrode %d',targi,ei))

