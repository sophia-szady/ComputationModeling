% Sophia Szady
% CS346 Computational Modeling
% HW 2 Exercise 1
% October 17, 2023

% Modeling a Lotka-Volterra Model between tuna (prey species Y), sharks
% (predator species P), and humans (H) that hunts both

% initial populations of the three species
init_Y_population = 100;
init_P_population = 15;
H_population = 1; % starting value at 1 so it doesn't impact the equations

% birth rates for the predator and prey
birth_rate_Y = 2;
birth_rate_P = 0.01;

% death rates for the predator and prey
death_rate_Y = 0.02;
death_rate_P = 1.06;

time_step = 0.001; % small time step for a more accurate model
num_steps = 24; % modeled in months so 24 months is 2 years

% creating arrays for the populations of the prey and predator over time 
Y_pop = zeros(num_steps/time_step,1); 
P_pop = zeros(num_steps/time_step,1);


% number of times the loop will be run
loop_iter = num_steps/time_step;

%running the simulation where Lotka-Volterra equations will be used to
%calculate the change in population in the predator-prey model
for j = 1:75:301
    %increasing the initial prey population from 1-200 by intervals of 50
    init_Y_population = j;
    Y_pop(1) = init_Y_population;
    P_pop(1) = init_P_population;
    for i = 2:loop_iter
       Y_change = ((birth_rate_Y*Y_pop(i-1)) - (death_rate_Y*P_pop(i-1)*Y_pop...
       (i-1) + death_rate_Y*Y_pop(i-1)*H_population))*time_step;
       P_change = ((birth_rate_P*Y_pop(i-1)*P_pop(i-1)*H_population)... 
       - (death_rate_P*P_pop(i-1)*H_population))*time_step;
       Y_pop(i) = Y_pop(i-1) + Y_change;
       P_pop(i) = P_pop(i-1) + P_change;
    end
    
    % plotting the original parameter Lotka-Volterra model
    %figure
    plot(0:time_step:num_steps-time_step,Y_pop, 0:time_step:num_steps-...
    time_step, P_pop) 
    title(sprintf("Lotka Volterra when initial prey population = %d", j))
    legend("Prey","Predator")
    xlabel("Time (in months)")
    ylabel("Population")
end

for j = 1:5:25
    %increasing the initial prey population from 1-200 by intervals of 50
    H_population = j;
    for i = 2:loop_iter
       Y_change = ((birth_rate_Y*Y_pop(i-1)) - (death_rate_Y*P_pop(i-1)*Y_pop...
       (i-1) + death_rate_Y*Y_pop(i-1)*H_population))*time_step;
       P_change = ((birth_rate_P*Y_pop(i-1)*P_pop(i-1)*H_population)... 
       - (death_rate_P*P_pop(i-1)*H_population))*time_step;
       Y_pop(i) = Y_pop(i-1) + Y_change;
       P_pop(i) = P_pop(i-1) + P_change;
    end
    
    % plotting the original parameter Lotka-Volterra model
    figure
    plot(0:time_step:num_steps-time_step,Y_pop, 0:time_step:num_steps-...
    time_step, P_pop) 
    title(sprintf("Lotka Volterra when human population = %d", j))
    legend("Prey","Predator")
    xlabel("Time (in months)")
    ylabel("Population")
end
