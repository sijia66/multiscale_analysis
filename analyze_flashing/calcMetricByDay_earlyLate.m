function [M, SD, X] = calcMetricByDay_earlyLate(dat, dayID, metric, N)


[M_early, SD_early, x_early] = calcMetricByDay(dat, dayID, metric, N, 'first');

[M_late, SD_late, x_late] = calcMetricByDay(dat, dayID, metric, N, 'last');

if x_early ~= x_late
    error('something is wrong...')
end

M  = [M_early(:)'; M_late(:)'];
SD = [SD_early(:)'; SD_late(:)'];
X  = x_early;
