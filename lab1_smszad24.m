% Sophia Szady
% Lab 1: Runge-Kutta 4
% 11/2/23

% Exercise 2: Implementing a Runge-Kutta 4 simulation of an undamped
% weighted spring system

%declaring initial parameters for the simulation

g = 9.81; % gravitational constant in m/s^2
k = 10; % spring constant in N/m
m = 0.2; % mass of the weight on the spring 
unweighted_length = 1; % length of the spring without a weight in m
init_displacement = 0.3; % initial displacement of the spring in m
weight_displacement = m*g / k; % amount the weight pulls down on the spring
init_length = unweighted_length + init_displacement + weight_displacement; 
% position at the start of the simualation
displacement = init_length - unweighted_length; %distance from the original
% position
time_step = 0.02; % time step in the simulation
simulation_length = 3; % length of the simulation
num_steps = simulation_length/time_step; % number of iterations in the 
% simulation 
weight = g * m ; % weight = gravitational constant * mass

% setting up an array for the acceleration, velocity, and position values 
% for each time step in the simulation
e_accelerations = zeros(num_steps,1); 
e_velocities = zeros(num_steps,1); 
e_positions = zeros(num_steps,1);

% Euler implementation to understand the system before implementing RK4

e_velocities(1) = 0; % position is not changing at the beginning
e_positions(1) = init_length; % where the spring starts

dvdt = @(ts, pos, vel) ((-k * (pos-unweighted_length))+weight)/m; 
% the rate of change of velocity (the acceleration)
dpdt = @(ts, pos, vel) vel;
% the rate of change of the position (the acceleration)

% As this is a second order system the acceleration is calculated to
% determine the velocity and the velocity determines the position, and the
% position from the previous time step helps calculate the acceleration
for i = 2:num_steps
    % calculating the velocity and position at given time step i 
    e_velocities(i) = e_velocities(i-1) + dvdt(time_step, e_positions(i-1), ...
        e_velocities(i-1)) * time_step;
    e_positions(i) = e_positions(i-1) + dpdt(time_step, e_positions(i-1), ...
        e_velocities(i-1)) *time_step;
end

% graphing the position with respect to time to see the oscillation of the
% spring
figure
plot(0:time_step:simulation_length-time_step,e_positions)
xlabel("time (in seconds)")
ylabel("position (in meters)")
title("The Oscillation of a Spring using Euler's Method")

% RK4 Simulation for the undamped weighted spring system

% initializing the arrays for the velocities and positions over the course
% of the simulation
RK4_velocities = zeros(num_steps,1);
RK4_positions = zeros(num_steps,1);

RK4_accelerations(1) = (-k*displacement+weight)/m; % using initial values
RK4_velocities(1) = 0; % position is not changing at the beginning
RK4_positions(1) = init_length; % where the spring starts


for i = 2:num_steps
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
    delta_velocity4 = dvdt(time_step*(i-1)+time_step, RK4_positions(i-1)+ ...
        delta_position3, RK4_velocities(i-1)+delta_velocity3)*time_step;
    delta_position4 = dpdt(time_step*(i-1)+time_step, RK4_positions(i-1)+ ...
        delta_position3, RK4_velocities(i-1)+delta_velocity3)*time_step;
    % a weighted average is calculated that places more weight on delta 2
    % and delta 3 to calculate the estimated change in velocity and
    % position
    delta_velocity = ((delta_velocity1 + (2*delta_velocity2) + ...
        (2*delta_velocity3) + delta_velocity4)/6);
    delta_position = ((delta_position1 + (2* delta_position2) + ...
        (2*delta_position3) + delta_position4)/6);
    % the estimated change is added to the previous value of velocity and
    % position
    RK4_velocities(i) = RK4_velocities(i-1) + delta_velocity;
    RK4_positions(i) = RK4_positions(i-1) + delta_position;
end

% plotting the position over time for the simulation
figure
plot(0:time_step:simulation_length-time_step, RK4_positions)
xlabel("time (in seconds)")
ylabel("position (in meters)")
title("The Oscillation of a Spring using RK4")


