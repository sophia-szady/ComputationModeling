% CS346 Smaller Assignment 2 Exercise 1
% Sophia Szady and MJ Pennington
% November 9th

% Implementing a model of a ball being thrown straight up from the side of
% a bridge using Euler's method

init_pos = 11; % initial position from the book in meters
init_vel = 15; % initial velocity from the book in meters/second
g = -9.81; % gravitational constant  
simulation_length = 4; % length of simulation in the textbook
time_step = 0.25; % length of time between each calculation
num_steps = simulation_length/time_step; % number of times position and 
% velocity are calculated 

% initializing arrays to store the values at each time step
velocities = zeros(num_steps,1);
positions = zeros(num_steps,1);

% setting the initial values to be the values from the textbook
velocities(1) = init_vel;
positions(1) = init_pos;

% anonymous functions to calculate the rate of change of velocity (dvdt)
% and the rate of change of position (dpdt)
dvdt = @(ts,p,v) g;
dpdt = @(ts,p,v) v;

% running the simulation for the specified length using Euler's method
for i = 2:num_steps+1
    % the new velocity and position values are calculated by determining
    % the rate of change, multiplying it by the time step and then adding
    % it to the previous value
    velocities(i) = velocities(i-1) + dvdt(time_step, positions(i-1), ...
        velocities(i-1)) * time_step;
    positions(i) = positions(i-1) + dpdt(time_step, positions(i-1), ...
        velocities(i-1)) *time_step;
end

% plotting the positions and velocities of the simulation
figure
plot(0:time_step:simulation_length,positions,0:time_step:simulation_length,velocities)
xlabel("time (in seconds)")
ylabel("position (in meters)")
title("Tossing a Ball off of a Bridge using Euler's Method")

