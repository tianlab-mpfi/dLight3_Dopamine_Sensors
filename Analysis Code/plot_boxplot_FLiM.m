%Calculates plot the boxplot of the FLIM experiment
clear all
close all

titles = {'agonist' 'antagonist'}

agonist = {'1' 2.49756618
'2'  2.563574737
'3' 2.562094873
'4' 2.558629533
'5' 2.533998181
'6' 2.530627033
'7'  2.543795718
'8' 2.550589722
'9' 2.560359007
'10' 2.542062987
};

antagonist = {'1' 2.304256674
'2'  2.325287911
'3' 2.373101825
'4' 2.331619241
'5' 2.335127966
'6' 2.324003059
'7'  2.334194504
'8' 2.336702043
'9' 2.332383889
'10' 2.328543912
};

combinedData = {agonist antagonist};

figure()

parrallel = [1 1 1 1];

for z = 1:2
    
    %% Load some sample data:
    figure()
    data = combinedData{z};
    
    measures = cell2mat(data(:,2));


    coordLineStyle = 'k.';
    boxplot(measures(1:size(data,1),1), 'Symbol', coordLineStyle);
    hold on
    scatter(1,measures(1:size(data,1),1), 'Color', 'k', 'Marker', '.')


    % if parrallel(z) == 1
    %     parallelcoords(measures(1:size(data,1),1:2), 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
    %       'Marker', '.', 'MarkerSize', 10, 'LineWidth',1);
    %         hold on
    % end
    ylim([2.29 2.58])

    title(titles{z})

end







