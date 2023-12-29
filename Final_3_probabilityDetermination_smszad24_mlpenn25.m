% Final Project File 3: Determining the simulation probabilities with 5
% randomly placed spores and Moore neighborhoods
% Sophia Szady and MJ Pennington
% 12/4/23

% We will be modeling the growth of fairy rings of mushrooms, a phenomena
% that describes the growth of mushrooms in circles

% This simulation will have 5 spores randomly placed in the field, and will
% be trying to determine probabilities for the model. A Moore neighborhood 
% and absorbing boundary conditions with a constant value of empty (0) will
% be used

% Potential cell values for cellular automata model of fairy rings that 
% represent the stage of the mushroom life cycle

empty = 0; % light green in the visualization, if any of the neighbors are
% young and probSpread determines if the cell stays empty of becomes a 
% young hyphae (young,2)
spore = 1; % black in the visualization, will either become a young hyphae
% (young,2) or stay a spore (1) depending on probSporeToHyphae
young = 2; % dark gray in the visualization, will become maturing (3) in 
% the next time step
maturing = 3; % light gray in visualization, will become mushrooms (4) or
% older depending on probMushroom
mushrooms = 4; % white in the visualization, will become decaying (6) in
% in next time step
older = 5; % light gray in the visualization, will become decaying (6) in 
%  the next time step
decaying = 6; % tan in the visualization, will become dead1 (7) in the next 
% time step
dead1 = 7; % brown in the visualization, will become dead2 (8) in the next 
% time step
dead2 = 8; % dark green in the visualization, will become empty (0) in the
% next time step

% dimensions of the field in the simulation r*c
r = 50;
c = 100;

% setting up the initial field and field with the boundary
initField = zeros(r,c);
initExtField = zeros(r+2,c+2);
%dField = zeros(x,y);

% initializing random coordinates for the spores
initSporeX1 = randi(r);
initSporeY1 = randi(c);
initSporeX2 = randi(r);
initSporeY2 = randi(c);
initSporeX3 = randi(r);
initSporeY3 = randi(c);
initSporeX4 = randi(r);
initSporeY4 = randi(c);
initSporeX5 = randi(r);
initSporeY5 = randi(c);

% setting the position of the initial spores
initField(initSporeX1,initSporeY1) = spore;
initField(initSporeX2,initSporeY2) = spore;
initField(initSporeX3,initSporeY3) = spore;
initField(initSporeX4,initSporeY4) = spore;
initField(initSporeX5,initSporeY5) = spore;
initExtField(2:r+1,2:c+1) = initField;

% amount of times the cells will be updated
numIterations = 40;

% cell arrays used to store each iterations of the field simulation,
% extFieldList includes enough cells for a boundary
fieldList{numIterations} = [];
extFieldList{numIterations} = [];

% making a cell array of the field and extended field, starting with the
% first time step, a new field and extended field will be added in the 
% simulation loop during each iteration
fieldList{1} = initField;
extFieldList{1} = initExtField;

% probabilities from one stage to another, all probabilites not included 
% are 1, meaning that they will always go to the next stage in the 
% mushroom life cycle
probSporeToHyphae = 0.5; % from spore (1) to young (2), if not it stays
% a spore
probMaturingToMushroom = 0.5; % from maturing (3) to mushroom (4), if not
% it becomes older
probSpread = 0.5; % from empty (0) to young (2), one of the neighbors needs
% to be young as well, if not it stays empty

% probabilities determined to demonstrate signficant changes in the
% simulation when the probabilities, probSporeToHyphae, 
% probMaturingToMushroom, or probSpread are changed one at a time
probsToTest=[0.1 .5 .9];

% an anonymous function to get the current value of the neighbors, using
% the extended field in order to account for the cells on the boundary.
% In order to get the right indices in the field, +1 is added to currentX
% and currentY to offset the extField
mooreNeighborVals = @(currentR, currentC, extField) ...
[extField(currentR+1+1,currentC+1) extField(currentR+1+1,currentC+1+1) ...
 extField(currentR+1,currentC+1+1) extField(currentR+1-1,currentC+1+1) ...
 extField(currentR+1-1,currentC+1) extField(currentR+1-1,currentC+1-1) ...
 extField(currentR+1,currentC+1-1) extField(currentR+1+1,currentC+1-1)];

% the simulation loop runs through all the cells and updates them based on 
% the progression of the mushroom life cycle, with some stages based on 
% probabilties, while others will always move to the next stage
for l=1:length(probsToTest) % going through each probability in the 
    % pre-selected probsToTest array for probSporeToHyphae
    probSporeToHyphae=probsToTest(l);
    for i = 2:numIterations % starting after the initial field setup
        field = fieldList{i-1}; % setting the field to the previous field
        extField = extFieldList{i-1}; % setting extField to the previous 
        % field
        for m = 1:r % going through each row in the field
            for n= 1:c % going through each column in the field
                current = field(m,n); % setting the current cell value
                neighborList = mooreNeighborVals(m,n, extField); % getting
                % the 8 neighbors of the current cell 
                if current == empty 
                    if ismember(young,neighborList) % if any of the
                    % neighbors are young
                        emptyRand = randi(10); % randomly calculating a
                        % value
                        if emptyRand < 10*probSpread % comparing the spread
                        % probability to the random number
                            field(m,n) = young; % if the value is greater 
                            % than the probability it becomes young, if not
                            % it stays empty
                        end
                    end
                elseif current == spore 
                    sporeRand = randi(10); % randomly calculating a value
                    if sporeRand < 10*probSporeToHyphae % comparing the 
                        % hyphae probability to the random number
                        field(m,n) = young; % if the value is greater than
                        % the probability it becomes young, if not it stays
                        % a spore
                    end
                elseif current == maturing 
                    maturingRand = randi(10); % randomly calculating a 
                    % value
                if maturingRand < 10*probMaturingToMushroom % comparing the 
                    % mushroom probability to the random number
                    field(m,n) = mushrooms; % if the value is greater than
                    % the probability it becomes mushrooms
                else
                    field(m,n) = older; % if the value is less than the
                    % probability it becomes older
                end
                elseif current == dead2
                    field(m,n) = empty; % if the mushroom is dead2 it 
                    % becomes empty 
                else
                    field(m,n) = current + 1; % if the current cell is not 
                    % one of the above special cases it moves on to the 
                    % next stage in the mushroom life cycle 
                end
            end
        end
    extField(2:r+1,2:c+1) = field; % setting the current field to be 
    % surrounded by a row of cells on each side to represent the boundary
    % conditions with the boundary cells being equal to empty (0)
    fieldList{i} = field; % saving the current fields in the list of all
    % the fields for visualization
    extFieldList{i} = extField; % saving the current extField in the list 
    % of all the fields for visualization
    end
    % a custom color map created to reflect the color values of each state 
    % of the mushroom life cycle, on the 0-1 rgb scale 
    lifeCycleMap = [ 0.702 0.98 0.702 % empty or 0 is light green  
    0 0 0 % spore or 1 is black
    0.349 0.349 0.349 % young or 2 is dark gray
    0.839 0.839 0.839 % maturing or 3 is light gray
    1 1 1 % mushrooms or 4 is white
    0.839 0.839 0.839 % older or 5 is light gray
    1 0.694 0.294 % decaying or 6 is tan
    0.478 0.318 0.102 % dead1 or 7 is brown
    0.129 0.349 0.133]; % dead2 or 8 is dark green 
    
    % visualizing each field in the simulation through a different frame
    for i=1:length(fieldList)
        fieldData = fieldList{i};
        colormap(lifeCycleMap); % setting the color map
        imagesc(fieldData); % visualizing the field
 
        % creating a color bar in order to see which color corresponds with
        % each color using all the states of the mushroom life cycle
        caxis([0,8]);
        lifeCycleColors=colorbar;
        lifeCycleColors.Ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5];    
        lifeCycleColors.TickLabels={'empty','spore','young','maturing',...
        'mushrooms','older','decaying','dead1','dead2'};
  
        hold;
    
        axis equal; axis tight; axis xy;
        probTitle=sprintf(...
        "Mushroom Ring Visualization with ProbSporeToHyphae set to %f",...
        probSporeToHyphae);
        title(probTitle);
        % wait to go on to next image
        fprintf('Waiting for any key to be pressed\n'); % key press 
        % activated frame change
        w = waitforbuttonpress;
    end       
end

% the simulation loop runs through all the cells and updates them based on 
% the progression of the mushroom life cycle, with some stages based on 
% probabilties, while others will always move to the next stage
for l=1:length(probsToTest) % going through each probability in the 
    % pre-selected probsToTest array for probSporeToHyphae
    probMaturingToMushroom=probsToTest(l);
    
    for i = 2:numIterations % starting after the initial field setup
        field = fieldList{i-1}; % setting the field to the previous field
        extField = extFieldList{i-1}; % setting extField to the previous 
        % field
        for m = 1:r % going through each row in the field
            for n= 1:c % going through each column in the field
                current = field(m,n); % setting the current cell value
                neighborList = mooreNeighborVals(m,n, extField); % getting
                % the 8 neighbors of the current cell 
                if current == empty 
                    if ismember(young,neighborList) % if any of the 
                    % neighbors are young
                        emptyRand = randi(10); % randomly calculating a
                        % value
                        if emptyRand < 10*probSpread % comparing the spread
                        % probability to the random number
                            field(m,n) = young; % if the value is greater 
                            % than the probability it becomes young, if not
                            % it stays empty
                        end
                    end
                elseif current == spore 
                    sporeRand = randi(10); % randomly calculating a value
                if sporeRand < 10*probSporeToHyphae % comparing the 
                % hyphae probability to the random number
                    field(m,n) = young; % if the value is greater than
                    % the probability it becomes young, if not it stays
                    % a spore
                end
                elseif current == maturing
                    maturingRand = randi(10); % randomly calculating a 
                    % value
                if maturingRand < 10*probMaturingToMushroom % comparing the 
                    % mushroom probability to the random number
                    field(m,n) = mushrooms; % if the value is greater than
                    % the probability it becomes mushrooms
                else
                    field(m,n) = older; % if the value is less than the
                    % probability it becomes older
                end
                elseif current == dead2
                    field(m,n) = empty; % if the mushroom is dead2 it 
                    % becomes empty 
                else
                    field(m,n) = current + 1; % if the current cell is not 
                    % one of the above special cases it moves on to the 
                    % next stage in the mushroom life cycle
                end
            end
        end
    extField(2:r+1,2:c+1) = field; % setting the current field to be 
    % surrounded by a row of cells on each side to represent the boundary
    % conditions with the boundary cells being equal to empty (0)
    fieldList{i} = field; % saving the current fields in the list of all
    % the fields for visualization
    extFieldList{i} = extField; % saving the current extField in the list 
    % of all the fields for visualization
    end
    
    % a custom color map created to reflect the color values of each state 
    % of the mushroom life cycle, on the 0-1 rgb scale 
    lifeCycleMap = [ 0.702 0.98 0.702 % empty or 0 is light green  
    0 0 0 % spore or 1 is black
    0.349 0.349 0.349 % young or 2 is dark gray
    0.839 0.839 0.839 % maturing or 3 is light gray
    1 1 1 % mushrooms or 4 is white
    0.839 0.839 0.839 % older or 5 is light gray
    1 0.694 0.294 % decaying or 6 is tan
    0.478 0.318 0.102 % dead1 or 7 is brown
    0.129 0.349 0.133]; % dead2 or 8 is dark green 
    
    % visualizing each field in the simulation through a different frame
    for i=1:length(fieldList)
        fieldData = fieldList{i};
        colormap(lifeCycleMap); % setting the color map
        imagesc(fieldData); % visualizing the field
 
        % creating a color bar in order to see which color corresponds with
        % each color using all the states of the mushroom life cycle
        caxis([0,8]);
        lifeCycleColors=colorbar;
        lifeCycleColors.Ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5];    
        lifeCycleColors.TickLabels={'empty','spore','young','maturing',...
        'mushrooms','older','decaying','dead1','dead2'};

        hold;
    
        axis equal; axis tight; axis xy;
        probTitle=sprintf(...
       "Mushroom Ring Visualization with ProbMaturingToMushrooms set to %f"...
        ,probMaturingToMushroom);
        title(probTitle);
        % wait to go on to next image
        fprintf('Waiting for any key to be pressed\n'); % key press 
        % activated frame change
        w = waitforbuttonpress;
    end   
end

% the simulation loop runs through all the cells and updates them based on 
% the progression of the mushroom life cycle, with some stages based on 
% probabilties, while others will always move to the next stage
for l=1:length(probsToTest) % going through each probability in the 
    % pre-selected probsToTest array for probSporeToHyphae 
    probSpread=probsToTest(l);
    disp(l)
    
    for i = 2:numIterations % starting after the initial field setup
        field = fieldList{i-1}; % setting current field to the previous 
        % field 
        extField = extFieldList{i-1}; % setting extField to the previous 
        % field
        for m = 1:r % going through each row in the field
            for n= 1:c % going through each column in the field
                current = field(m,n); % setting the current cell value
                neighborList = mooreNeighborVals(m,n, extField); % getting
                % the 8 neighbors of the current cell 
                if current == empty 
                    if ismember(young,neighborList) % if any of the 
                    % neighbors are young
                        emptyRand = randi(10); % randomly calculating a
                        % value
                        if emptyRand < 10*probSpread % comparing the spread
                        % probability to the random number
                            field(m,n) = young; % if the value is greater 
                            % than the probability it becomes young, if not
                            % it stays empty
                        end
                    end
                elseif current == spore 
                    sporeRand = randi(10); % randomly calculating a value
                if sporeRand < 10*probSporeToHyphae % comparing the 
                % hyphae probability to the random number
                    field(m,n) = young; % if the value is greater than
                    % the probability it becomes young, if not it stays
                    % a spore
                end
                elseif current == maturing
                    maturingRand = randi(10); % randomly calculating a 
                    % value
                if maturingRand < 10*probMaturingToMushroom % comparing the 
                    % mushroom probability to the random number
                    field(m,n) = mushrooms; % if the value is greater than
                    % the probability it becomes mushrooms
                else
                    field(m,n) = older; % if the value is less than the
                    % probability it becomes older
                end
                elseif current == dead2
                    field(m,n) = empty; % if the mushroom is dead2 it 
                    % becomes empty 
                else
                    field(m,n) = current + 1; % if the current cell is not 
                    % one of the above special cases it moves on to the 
                    % next stage in the mushroom life cycle
                end
            end
        end
    extField(2:r+1,2:c+1) = field; % setting the current field to be 
    % surrounded by a row of cells on each side to represent the boundary
    % conditions with the boundary cells being equal to empty (0)
    fieldList{i} = field; % saving the current fields in the list of all
    % the fields for visualization
    extFieldList{i} = extField; % saving the current extField in the list 
    % of all the fields for visualization
    end
    
    % a custom color map created to reflect the color values of each state of
    % the mushroom life cycle, on the 0-1 rgb scale 
    lifeCycleMap = [ 0.702 0.98 0.702 % empty or 0 is light green  
    0 0 0 % spore or 1 is black
    0.349 0.349 0.349 % young or 2 is dark gray
    0.839 0.839 0.839 % maturing or 3 is light gray
    1 1 1 % mushrooms or 4 is white
    0.839 0.839 0.839 % older or 5 is light gray
    1 0.694 0.294 % decaying or 6 is tan
    0.478 0.318 0.102 % dead1 or 7 is brown
    0.129 0.349 0.133]; % dead2 or 8 is dark green 
    
    % visualizing each field in the simulation through a different frame
    for i=1:length(fieldList)
        fieldData = fieldList{i};
        colormap(lifeCycleMap); % setting the color map
        imagesc(fieldData); % visualizing the field
 
        % creating a color bar in order to see which color corresponds with
        % each color using all the states of the mushroom life cycle
        caxis([0,8]);
        lifeCycleColors=colorbar;
        lifeCycleColors.Ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5];    
        lifeCycleColors.TickLabels={'empty','spore','young','maturing',...
        'mushrooms','older','decaying','dead1','dead2'};
  
        hold;
    
        axis equal; axis tight; axis xy;
        probTitle=sprintf(...
        "Mushroom Ring Visualization with Prob Spread set to %f", ...
        probSpread );
        title(probTitle);
        % wait to go on to next image
        fprintf('Waiting for any key to be pressed\n'); % key press 
        % activated frame change
        w = waitforbuttonpress;
    
  
    end
 
        
end

        
        
        