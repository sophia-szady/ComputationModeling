% Sophia Szady
% CS346 Computational Modeling
% HW 1 Exercise 1
% September 25, 2023

% quantity of radium-226 (in percent)
Q0 = 1; 
% amount of time the simulation (in years)
totalTime = 10000;
% time step (in years)
step = 0.5;
% rate of decay of radium-226 (in percent)
decayRate = 0.0427869;

% number of times the simulation will be run
numIter = totalTime / step;

over60 = 0; % 0 indicates the percent left is over 60
allQ0 = ones(numIter+1);

%running the simulation for the given # of iterations
for i = 1:numIter
    allQ0(i) = Q0;
    growthRate = decayRate * Q0; % calculating the rate of decay
    popChange = growthRate * step; % calculating the amount of decay
    Q0 = Q0 - popChange; % subtracting the decay from the amount left
    t = i * step; % calculating the current step
    if t == 500 % percent left after 500 years 
        disp(Q0);
    end
    if t == 5000 % percent left after 5000 years
        disp(Q0);
    end
    if Q0 <= .6 && over60 == 0 % the percent remaining drops under 60% 
                               % for the first time
        disp(t);
        over60 = 1; %making it so the loop isn't entered again
    end
end

figure % creating a new plot for the function
plot((0:step:totalTime), allQ0) % graphing Q0 at each time step
title('Decay of Radium-226 over 20,000 years') % titling the graph
xlabel('Time (in years)')
ylabel('Percent remaining of Radium-226')