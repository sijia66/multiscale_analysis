function plotTF_dB(s,freqRaw,freqRange,trLen)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s_length = size(s,1);

s_timePoints = trLen * 0.5 +(1:s_length)*trLen * 0.1;
s_step = 10;



freqRangeIndex = cal_index_freq(freqRaw ,freqRange(1),freqRange(2));
freq = freqRaw(freqRangeIndex);

freq_length = length(freq);
freq_step  = round(freq_length / 32);


imagesc((s(:,1:freq_length)'))
colorbar
%caxis([-2 2])
set(gca,'YDir','normal')
set(gca, 'XTick', 5:s_step:s_length, 'XTickLabel', s_timePoints(5:s_step:s_length))
set(gca, 'YTick', [1:freq_step:freq_length], 'YTickLabel', round(freq(1:freq_step:freq_length))) ;% 20 ticks
ylabel('Frequency (Hz)')
xlabel('Time (s)')

end

