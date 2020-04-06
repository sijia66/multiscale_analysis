%%
% we just wanna look at time domain
ecog_ei = 1;
sc32_ei = 1;

plot(squeeze(mean(dataTarget.target1(:,ecog_ei,:))))

hold on
plot(squeeze(mean(dataTarget_SC32.target1(:,ecog_ei,:))))
