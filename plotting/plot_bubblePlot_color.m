



freqRange_display = [5 13];
freqRange = find(freq <= freqRange_display(2) & freq >= freqRange_display(1));
bubbleValues = mean(mi_x_y(:,freqRange),2);

figure

magFac = 200;
labelElectrode = 1;
textShift = 0.1;
titleTxt = 'Mutual Information With Target Reach Direction';


bubbleValues_inc_indices = find(bubbleValues > 0 );
bubbleValues_dec_indices = find(bubbleValues < 0);
bubbleValues_zero_indices = find(bubbleValues == 0);

bubbleValues = bubbleValues / max(bubbleValues);

hold on
%plot positive values
for bi = 1:length(bubbleValues_inc_indices)
    
scatter(P_chamber(bubbleValues_inc_indices(bi),1),P_chamber(bubbleValues_inc_indices(bi),2)...
    ,magFac, ...
    'MarkerFaceColor',[bubbleValues(bubbleValues_inc_indices(bi)) 0 0])
axis equal

end

%plot negative values
scatter(P_chamber(bubbleValues_dec_indices,1),P_chamber(bubbleValues_dec_indices,2)...
    ,abs(bubbleValues(bubbleValues_dec_indices)) * magFac,...
    'MarkerFaceColor','b','MarkerEdgeColor','none')

scatter(P_chamber(bubbleValues_zero_indices,1),P_chamber(bubbleValues_zero_indices,2)...
    ,magFac  / 2, 'MarkerFaceColor','g')
hold off

title(titleTxt)
xlabel('mm')
ylabel('mm')

if labelElectrode
    
    hold on
    if drive_idx == 1 %ECOG
        for i = 1:size(bubbleValues,1)
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
