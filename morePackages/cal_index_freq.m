function freqRange = cal_index_freq(specDatFreq,lowFreq,highFreq)
% this function returns inclusive frequency range index

tempRange = find(specDatFreq >= lowFreq);
lowerBound = tempRange(1);

tempRange = find(specDatFreq <= highFreq);
highBound = tempRange(end);

freqRange = lowerBound:highBound;

end

