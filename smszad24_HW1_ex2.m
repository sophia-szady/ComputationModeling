% Sophia Szady
% CS346 Computational Modeling
% HW 1 Exercise 2
% September 25, 2023

n = 40000; % number of points generated
r = 1; % radius of the unit circle 

x = rand([1,n]); % generating n number of random x values between 0 and 1
y = rand([1,n]); % generating n number of random y values between 0 and 1
z = x.^2 + y.^2; % calculating whether the point is inside the unit
% circle using the formula x^2 + y^2 = r for the area of the circle
inside = sum(z<r)/n; % if z is larger than r it means that the point is 
% outside the unit circle, summing and dividing by the number of points
% generated finds the percent of points inside the circle
disp(inside*4) % multiplying inside by 4 estimates pi as our calculations
% has only considered one coordinate

figure % creating a new figure to represent points inside and outside the
% unit circle 
plot(x(z<=r),y(z<=r), 'r.', x(z>r),y(z>r), 'b.') % plotting all of the simulated points by whether they 
% are inside the circle or not
legend("inside the circle","outside the circle"); % if the points is inside
title("Estimating Pi Using a Monte Carlo Simulation")
% the circle it is red and if it is outside it is blue 
max = 0; % setting a lower threshold that all estimates will be higher than
min = 5; % setting a higher threshold that all estimates will be higher than
avg = 0; % setting a value for the average pi estimate
numSim = 1000; % number of times pi will be estimated

%loop to estimate pi numSim times
for i = 1:numSim 
    % repeating the same calculations as lines 9-13
    x = rand([1,n]);
    y = rand([1,n]);
    z = x.^2 + y.^2;
    inside = sum(z<r)/n;
    avg = avg + (inside*4); % summing all the pi estimations
    if (inside*4) > max % if the current pi estimate is the biggest
        max = (inside*4);
    end
    if (inside*4) < min % if the current pi estimate is the smallest
        min = (inside*4);
    end
end
% dividing the sum of pi estimates by the number of estimates
avg = avg/numSim; 
fprintf("The minimum is: %f\n the max is: %f\n the average is: %f\n", ... 
    min, max, avg);