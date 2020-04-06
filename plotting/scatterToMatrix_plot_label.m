%scatterToMatrix
function scatterToMatrix_plot_label(fig, map, mapInfo)


axis xy
axis equal
scatterToMatrix_plot_step = 4;
%       only label x axis
[x_num,y_num] = size(map);
xticks = 1:scatterToMatrix_plot_step:x_num;
xlabels = ((0:scatterToMatrix_plot_step:x_num-1)) * mapInfo.x_pitch + mapInfo.x_zero;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xlabel('mm')


yticks = 1:scatterToMatrix_plot_step:y_num;
ylabels = ((0:scatterToMatrix_plot_step:x_num-1)) * mapInfo.y_pitch + mapInfo.y_zero;
set(gca, 'YTick', xticks, 'YTickLabel', xlabels);
ylabel('mm')
colorbar

set(gca, 'fontsize', 16)