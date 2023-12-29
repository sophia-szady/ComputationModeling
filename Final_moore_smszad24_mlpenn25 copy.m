% Final Project
% Sophia Szady and MJ Pennington
% 12/4/23

% We will be modeling the growth of fairy rings of mushrooms, a phenomena
% that describes the growth of mushrooms in circles

% This simulation will have a single spore randomly placed in the field,
% and uses a von Neumann neighborhood and absorbing boundary conditions 
% with a constant value of empty (0).

% Cell values for cellular automata model of fairy rings that represent ...
% the stage of the mushroom life cycle

empty = 0; % light green in the visualization, probYoung and if the ...
% neighbors are young determines if the cell stays empty of becomes a young
% hypahe
spore = 1; % black in the visualization, will either become a young hyphae
% or stay a spore depending on probSporeToHyphae
young = 2; % dark gray in the visualization, will become maturing in the
% next time step
maturing = 3; % light gray in visualization, will become mushrooms or ...
% older depending on probMushroom
mushrooms = 4; % white in the visualization, will become decaying in the
% next time step
older = 5; % light gray in the visualization, will become decaying in the 
% next time step
decaying = 6; % tan in the visualization, will become dead1 in the next 
% time step
dead1 = 7; % brown in the visualization, will become dead2 in the next 
% time step
dead2 = 8; % dark green in the visualization, will become empty in the next
% time step

% dimensions of the simulation r*c
r = 50;
c = 100;

% setting up the initial field and extField with the boundary condition
% cells
initField = zeros(r,c);
initExtField = zeros(r+2,c+2);

% initializing a random coordinate for the first spore
initSporeX1 = randi(r);
initSporeY1 = randi(c);

% setting the position of the initial spore
initField(initSporeX1,initSporeY1) = spore;
initExtField(2:r+1,2:c+1) = initField;

% amount of times the cells will be updated
numIterations = 100;

% cell arrays used to store each iterations of the field simulation,
% extFieldList includes enough cells for boundar
fieldList{numIterations} = [];
extFieldList{numIterations} = [];

% making a cell array of the field and extended field, starting with the
% first time step, a new field and extended field will be added in the 
% simulation loop during each iteration
fieldList{1} = initField;
extFieldList{1} = initExtField;

% probabilities from one stage to another, all not included are 1
probSporeToHyphae = 0.8;
probMaturingToMushroom = 0.7;
probSpread = 0.9;

% an anonymous function to get the current value of all the neighbors
vonNeumannNeighborVals = @(currentX, currentY, extField) ...
    [extField(currentX+1,currentY+1+1) extField(currentX+1+1, currentY+1)...
    extField(currentX+1-1, currentY+1) extField(currentX+1, currentY+1-1)];

% the simulation loop runs through all the cells and updates them based on 
% the progression of the mushroom life cycle, with some stages based on 
% probabilties, while others will move to the next stage
for i = 2:numIterations % starting after the initial field setup
    field = fieldList{i-1}; % setting current field to the previous field 
    extField = extFieldList{i-1}; % setting extField to the previous field
    for m = 1:r % going through each row in the field
        for n= 1:c % going through each column in the field
            current = field(m,n); % setting the current cell value
            neighborList = vonNeumannNeighborVals(m,n, extField); % getting
            % the 4 neighbors of the current cell 
            if current == empty 
                if ismember(young,neighborList) % if any of the neighbors
                    % are young
                    emptyRand = randi(10); % randomly calculating a value
                    if emptyRand < 10*probSpread % comparing the spread
                        % probability to the random number
                        field(m,n) = young; % if the value is greater than
                        % the probability it becomes young, if not it stays
                        % young
                    end
                end
            elseif current == spore 
                sporeRand = randi(10); % randomly calculating a value
                if sporeRand < 10*probSporeToHyphae % comparing the hyphae
                        % probability to the random number
                    field(m,n) = young; % if the value is greater than
                    % the probability it becomes young, if not it stays
                    % young
                end
            elseif current == maturing 
                maturingRand = randi(10); % randomly calculating a value
                if maturingRand < 10*probMaturingToMushroom % comparing the 
                    % mushroom probability to the random number
                    field(m,n) = mushrooms; % if the value is greater than
                    % the probability it becomes mushrooms
                else
                    field(m,n) = older; % if the value is less than the
                    % probability it becomes older
                end
            elseif current == dead2
                field(m,n) = empty; % if the mushroom is dead2 it becomes
                % empty 
            else
                field(m,n) = current + 1; % if the current cell is not one
                % of the above special cases it moves on to the next stage
                % in the mushroom life cycle 
            end
        end
    end
    extField(2:r+1,2:c+1) = field; % setting the current field to be 
    % surrounded by a row of cells on each side to represent the boundary
    % conditions
    fieldList{i} = field; % saving the current fields in the list of all
    % the fields
    extFieldList{i} = extField; % saving the current extField in the list 
    % of all the fields
end
% a custom color map created to reflect the color values of each state of 
% the mushroom life cycle, on the 0-1 rgb scale 
map = [ 0.702 0.98 0.702 % empty or 0 is light green  
    0 0 0 % spore or 1 is black
    0.349 0.349 0.349 % young or 2 is dark gray
    0.839 0.839 0.839 % maturing or 3 is light gray
    1 1 1 % mushrooms or 4 is white
    0.839 0.839 0.839 % older or 5 is light gray
    1 0.694 0.294 % decaying or 6 is tan
    0.478 0.318 0.102 % dead1 or 7 is brown
    0.129 0.349 0.133]; % dead2 or 8 is dark green 
% creating another color map in order to just visualize the mushrooms, 
% where all values are light green except for mushrooms or 4 which is white
mushroommap= [ 0.702 0.98 0.702
               0 0 0
               0.702 0.98 0.702
               0.702 0.98 0.702
               1 1 1
               0.702 0.98 0.702
               0.702 0.98 0.702
               0.702 0.98 0.702
               0.702 0.98 0.702];
% visualizing each field in the simulation through a different frame with
% all colors visible for each value
for i=1:length(fieldList)
    data = fieldList{i}; 
    colormap(map); % setting the color map
    % plot the image
    imagesc(data); % visualizing the field
    title(sprintf('Frame: %d', i)); % setting the title to match the field
    % iteration
 
    % creating a color bar in order to see which color corresponds with
    % each color using all the states of the mushroom life cycle
    caxis([0,9]) 
    c=colorbar
    c.Ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]    
    c.TickLabels={'empty','spore','young','maturing','mushrooms','older', ...
        'decaying','dead1','dead2'}
  
    hold;
    
    axis equal; axis tight; axis xy;
    % wait to go on to next image
    fprintf('Waiting for any key to be pressed\n'); % key press activated 
    % frame change
    w = waitforbuttonpress;
end
% visualizing just the mushroom in order to see the rings distinctly
for i=1:length(fieldList)
    data = fieldList{i};
    colormap(mushroommap); % setting the colors to be all light green
    % besides white to see the mushrooms
    % plot the image
    imagesc(data);
    title(sprintf('Visualization of Mushrooms at Frame: %d', i));
    hold;
    
    axis equal; axis tight; axis xy;
    % wait to go on to next image
    fprintf('Waiting for any key to be pressed\n'); % key press activated
    % frame change
    w = waitforbuttonpress;
end
        
        
        