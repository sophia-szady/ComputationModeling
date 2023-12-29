% Sophia Szady
% CS346 Computational Modeling
% HW 1 Exercise 3
% September 27, 2023

W = 20; % the number of random walks in the simulation
numSteps = 20; % the max number of steps-1 because initial placement counts
B = 5; % the boundary distance, if a walker reaches this value whether 
% positive or negative, the walker will be frozen here for the rest of the
% simulation
walkers = zeros(numSteps,W); % an array of all the walkers and
% their walks
firstUnfrozen = 'False'; % haven't created the unfrozen walkers variable yet
unfrozenWalkers = zeros(numSteps,1); % used to store the walks of unfrozen walkers

for walkNum = 1:W % looping through each walker
    walkers(1,walkNum) = 0; % setting the initial position to 0
    for stepNum = 2:numSteps % going through each walk 
        %step 1 is the initial position so j starts at 2
        if walkers(stepNum-1,walkNum) > B %if the previous step reached the 
            %upper threshold
            walkers(stepNum,walkNum) = walkers(stepNum-1,walkNum); % if the previous value reached the threshold
            % all of the rest of the values are the threshold value too
        elseif walkers(stepNum-1,walkNum) < -B % checking the same as the if statement 
            % but on the lower threshold
            walkers(stepNum,walkNum) = walkers(stepNum-1,walkNum); % same logic as line 21
        else
            walkers(stepNum,walkNum) = walkers(stepNum-1,walkNum) + randn; % calculating the 
            % step next step for the walk
        end
    end
    %disp(sum((walkers(:,walkNum)>-B)))
    if sum((walkers(:,walkNum)<B)) == numSteps && sum((walkers(:,walkNum)>-B)) == numSteps 
        % checking if the walker ever freezes, if not enter loop 
        if walkNum == false % first time entering the loop
            unfrozenWalkers = walkers(:,walkNum); % declaring the first 
            % unfrozen walkers
            firstUnfrozen = 'True'; % the unfrozen walkers variable was declared
        end
        disp(walkNum)
        unfrozenWalkers = cat(2,unfrozenWalkers,walkers(:,walkNum));
        % adding the each following unfrozen walker
    end
end
numFrozen = sum(sum(walkers>B) > 0) + sum(sum(walkers<-B) > 0); % number of walkers that reach B or -B
avgFrozenSteps = sum(sum(walkers>B)) + sum(sum(walkers<-B)); % total sum of all frozen walkers
% steps it takes to reach the threshold B
figure % creating a figure for the walker plot
plot(0:1:numSteps-1,walkers) % graphing the walk
xlabel("timestep") 
ylabel("position")
text(0.5,5,"The Walks") % showing what the lines represent


avgFrozenSteps = avgFrozenSteps/numFrozen; % calculating the average number
% of steps it takes to reach B or -B only for those that do
fprintf("The average number of steps when frozen walkers collide: %f\n", ...
    avgFrozenSteps)
fprintf("The average number of frozen walkers: %d\n", ...
    numFrozen)
figure % creating a new figure for only walkers that walked all 20 steps
plot(0:1:numSteps-1,unfrozenWalkers) % plotting only walkers that didn't reach B or -B
xlabel("timestep")
ylabel("position")
text(0.5,4,"The Walks") % see line 49
