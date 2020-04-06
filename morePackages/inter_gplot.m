function  inter_gplot(A,graphParam)
%UNTITLED Summary of this function goes here
%   A is a pseudo Adjacency matrix and is logical size N1 * N2
%   
%   Note that the labeling is only row/column number and the actual E 
%   number

Nelectrodes1 = size(A,1);
Nelectrodes2 = size(A,2);

%check if the length of the coordinate vectors match the matrix dimensions

if size(graphParam.P_chamber1,1) ~= Nelectrodes1
    error('First dimension doesn t match coordinate length')
end 
if size(graphParam.P_chamber2,1) ~= Nelectrodes2
    error('Second dimension doesn t match coordinate length')
end 


axis([-6 6 -6 6])
axis square

hold on
for xi = 1:Nelectrodes1
    for yi = 1:Nelectrodes2
        if A(xi,yi)
            
            %draw the two end points
            scatter(graphParam.P_chamber1(xi,1),graphParam.P_chamber1(xi,2),[],'y','filled','o');
            scatter(graphParam.P_chamber2(yi,1),graphParam.P_chamber2(yi,2),50,'b','filled','o');
            
            %label the two points
            text(graphParam.P_chamber1(xi,1)+ graphParam.textShift, ...
                 graphParam.P_chamber1(xi,2), ...
                 num2str(xi));
             
             
            text(graphParam.P_chamber2(yi,1)+graphParam.textShift, ...
                 graphParam.P_chamber2(yi,2), ...
                 num2str(yi));
            
            %draw a line in between
            line([graphParam.P_chamber1(xi,1) graphParam.P_chamber2(yi,1)],... 
                 [graphParam.P_chamber1(xi,2) graphParam.P_chamber2(yi,2)],....
                 'color','k');
             
            
        end
    end
end

txt = sprintf('Cutoff Threshold %0.2f (%d th Percentile)',graphParam.threshold,graphParam.thresholdPercentile);
text(-3,-8,txt,'FontSize',14)

txt = sprintf('Averaged Between %0.2f and %0.2f Hz',graphParam.freqRange(1),graphParam.freqRange(2));
text(-3,-8.8,txt,'FontSize',14)

xlabel('mm')
ylabel('mm')

hold off

end

