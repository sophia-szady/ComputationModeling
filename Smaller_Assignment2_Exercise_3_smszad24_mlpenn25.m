% CS346 Smaller Assignment 2 Exercise 3
% Sophia Szady and MJ Pennington
% November 9th

% Implementing a model of a ball being thrown straight up from the side of
% a bridge using RK4 and a new acceleration equation

init_pos = 11; % initial position from the book in meters
init_vel = 15; % initial velocity from the book in meters/second
g = -9.81; % gravitational constant
simulation_length = 4; % length of simulation in the textbook
time_step = 0.25; % length of time between each calculation
num_steps = simulation_length/time_step; % number of times position and 
% velocity are calculated 

% initializing arrays to store the values at each time step
RK4_velocities = zeros(num_steps,1);
RK4_positions = zeros(num_steps,1);

% setting the initial values to be the values from the textbook
RK4_velocities(1) = init_vel;
RK4_positions(1) = init_pos;

% anonymous functions to calculate the rate of change of velocity (dvdt)
% and the rate of change of position (dpdt)
dvdt = @(ts,p,v) g+0.01*(v+p)+0.3*ts^2;
dpdt = @(ts,p,v) v;

% running the simulation for the specified length using RK4
for i = 2:num_steps+1
    % Delta 1 uses euler's method to calculate the change in position and 
    % velocity at the given time step
    delta_velocity1 = dvdt(time_step*(i-1), RK4_positions(i-1), ...
        RK4_velocities(i-1))* time_step; 
    delta_position1 = dpdt(time_step*(i-1), RK4_positions(i-1), ... 
        RK4_velocities(i-1))* time_step; 
    % Delta 2 estimates halfway between the previous time step and the
    % current time step using delta 1
    delta_velocity2 = dvdt(time_step*(i-1)+time_step*0.5, RK4_positions(i-1)+ ...
        delta_position1*0.5, RK4_velocities(i-1)+delta_velocity1*0.5)*time_step;
    delta_position2 = dpdt(time_step*(i-1)+time_step*0.5, RK4_positions(i-1)+ ...
        delta_position1*0.5, RK4_velocities(i-1)+delta_velocity1*0.5)*time_step;
    % Delta 3 estimates halfway between the previous time step and the
    % current time step using delta 2
    delta_velocity3 = dvdt(time_step*(i-1)+time_step*0.5, RK4_positions(i-1)+ ...
        delta_position2*0.5, RK4_velocities(i-1)+delta_velocity2*0.5)*time_step;
    delta_position3 = dpdt(time_step*(i-1)+time_step*0.5, RK4_positions(i-1)+ ...
        delta_position2*0.5, RK4_velocities(i-1)+delta_velocity2*0.5)*time_step;
    % Delta 4 estimates the change in position and velocity using delta 3 
    delta_velocity4 = dvdt(time_step*(i-1), RK4_positions(i-1)+ ...
        delta_position3, RK4_velocities(i-1)+delta_velocity3)*time_step;
    delta_position4 = dpdt(time_step*(i-1), RK4_positions(i-1)+ ...
        delta_position3, RK4_velocities(i-1)+delta_velocity3)*time_step;
    % a weighted average is calculated that places more weight on delta 2
    % and delta 3 to calculate the estimated change in velocity and
    % position
    delta_velocity = ((delta_velocity1 + (2*delta_velocity2) + ...
        (2*delta_velocity3) + delta_velocity4)/6);
    delta_position = ((delta_position1 + (2* delta_position2) + ...
        (2*delta_position3) + delta_position4)/6);
    disp(delta_position)
    % the estimated change is added to the previous value of velocity and
    % position
    RK4_velocities(i) = RK4_velocities(i-1) + delta_velocity;
    RK4_positions(i) = RK4_positions(i-1) + delta_position;
end

% plotting the positions and velocities of the simulation
figure
plot(0:time_step:simulation_length, RK4_positions,0:time_step:simulation_length,RK4_velocities)
xlabel("time (in seconds)")
ylabel("position (in meters)")
title("Simulating Throwing a Ball off a Bridge")

