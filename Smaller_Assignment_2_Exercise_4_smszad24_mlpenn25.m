% CS346 Smaller Assignment 2 Exercise 4
% Sophia Szady and MJ Pennington
% November 9th

% Implementing a model of a ball dropped down from the side of a bridge 
% using RK4

mass = 0.5; % initial mass of the ball from the textbook
g = -9.81; % gravitational constant
radius = 0.05; % initial radius of the ball from the textbook
weight = mass*g; % using mass and the gravitational constant to calculate 
%weight
p = pi; % using the built in matlab constant of pi
area = p * radius^2; % area of the ball
init_pos = 400; % initial position from the book in meters
init_vel = 0; % initial velocity from the book in meters/second
init_speed = abs(init_vel); % calculating the speed at the beginning of the
% simulation
init_air_friction = -0.65*area*init_vel*init_speed; % calculating the air
% friction at the beginning of the simulation
simulation_length = 15; % length of simulation in the textbook
time_step = 0.01; % the time between each calculation
num_steps = simulation_length/time_step; % number of times position, speed, 
% and velocity are calculated 

% initializing arrays to store the values at each time step
RK4_velocities = zeros(num_steps,1);
RK4_speeds = zeros(num_steps,1);
RK4_positions = zeros(num_steps,1);

% setting the initial values to be the values from the textbook
RK4_velocities(1) = init_vel;
RK4_positions(1) = init_pos;
RK4_speeds(1) = init_speed;

% anonymous functions to calculate the rate of change of velocity (dvdt)
% and the rate of change of position (dpdt)
dvdt = @(ts,p,v) (weight+ (-0.65*area*v*abs(v)))/mass;
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
    % the estimated change is added to the previous value of velocity and
    % position, the speed is the absolute value of velocity
    RK4_velocities(i) = RK4_velocities(i-1) + delta_velocity;
    RK4_speeds(i) = abs(RK4_velocities(i));
    RK4_positions(i) = RK4_positions(i-1) + delta_position;
end

% plotting the positions and velocities of the simulation
figure
plot(0:time_step:simulation_length, RK4_positions,0:time_step:simulation_length,RK4_speeds)
xlabel("time (in seconds)")
ylabel("position (in meters)")
title("Simulating Dropping a Ball off a Bridge")

