% HW 3 Exercise 4
% Sophia Szady and MJ Pennington
% 11/12/23

% Declaring variables for the simulation of the model
capacitance = 0.1; % ability to hold charge
I_ext = 15; % external current
init_n = 0.317; % potassium activation gating variable
init_m = 0.05; % sodium activation gating variable
init_h = 0.6; % sodium inactivation gating variable
gk = 36; % maximum potassium (K) conductance
gna = 120; % maximum sodium (Na) conductance
gl = 0.3; % maximum leakage conductance

init_V = -65; % voltage in mV
Vl = -54.4; % displacement from the equilibrium potential for leakage
Vk = -77; % displacement from the equilibrium potential for potassium
Vna = 50; % displacement from the equilibrium potential for sodium

% gate thresholds in mV
na_open = -55; % the sodium gate opens at this value
na_close = 50; % the sodium gate closes at this value
k_open = 50; % the potassium gate opens at this value
k_close = -75; % the potassium gate closes at this value

% initial concentrations in mM/L
Ki = 150; % initial potassium concentration inside the cell
Ko = 5.5; % initial potassium concentration outside the cell
Nai = 15; % initial sodium concentration inside the cell
Nao = 150; % initial sodiium concentration outside the cell

% opening rate constant for potassium activation
alpha_n = @(V,n,m,h) (0.01*(V+55))/(1-(exp(-(V+55)/10)));
% opening rate constant for sodium activation
alpha_m = @(V,n,m,h) (0.1*(V+40))/(1-(exp(-(V+40)/10)));
% opening rate constant for sodium inactivation
alpha_h = @(V,n,m,h) 0.07*exp(-((V+65)/20));
% closing rate constant for potassium activation
beta_n = @(V,n,m,h) 0.125*exp(-((V+65)/80));
% closing rate constant for sodium activation
beta_m = @(V,n,m,h) 4*exp(-((V+65)/18));
% closing rate constant for sodium inactivation
beta_h = @(V,n,m,h) 1/((exp(-(V+35)/10)+1));

time_step = 0.001; % in ms
simulation_length = 3; % length of the simulation in ms
ext_time = 0.5; % time when external current is applied
ext_length = 0.5; % time of external current being applied
total_steps = simulation_length/time_step; % 
null_steps = ext_time/time_step; % time step when external current is applied
step_off = (null_steps) + (ext_length/time_step); % the step number when the 
% applied current stops

% binary constants for potassium (k) and sodium (na) that coincide with the
% voltage gating, starting vals are 0 because the gates are closed at the
% beginning of the simulation
open_val_k = 0;
open_val_na = 0;

% anonymous functions that represent the potassium (Ik), sodium (Ina), and
% leakage (Il) currents 
% open_k and open_na are used to indicate when the gate is open, as the 
% whole term will be multiplied by 0 when the gate is closed and 1 when the
% gate is open
Ik = @(V,n,m,h,open_k) open_k*gk*(n^4)*(V-Vk); 
Ina = @(V,n,m,h,open_na) open_na*gna*(m^3)*h*(V-Vna);
Il = @(V,n,m,h) gl*(V-Vl);
Ip = Il(init_V,init_n,init_m,init_h);

% the Na-K pump moves 3 sodium ions out of the cell for every 2 potassium 
% ions that go into the cell
Na_flow = Ip * 3;
K_flow = Ip * 2;

% arrays for the action potential (V), potassium activation gating variable
% (n), sodium activation gating variable (m), and the sodium inactivation
% gating variable (h), for each time step in the simulation
V_vals = zeros(total_steps,1);
n_vals = zeros(total_steps,1);
m_vals = zeros(total_steps,1);
h_vals = zeros(total_steps,1);

% a loop setting all the values before the external current is applied to
% be the initial values
for i = 1:null_steps
    V_vals(i) = init_V;
    n_vals(i) = init_n;
    m_vals(i) = init_m;
    h_vals(i) = init_h;
end

% anonymous functions for the rates of change of the action potential (V),
% potassium activation gating variable (n), sodium activation gating 
% variable (m), and the sodium inactivation gating variable (h) 

% These functions use other anonymous functions, including the ones for 
% current, opening rate constant, and closing rate constant
% current = (I_ext-Ik(V,n,m,h,open_k)-Ina(V,n,m,h,open_na)-(Il(V,n,m,h)-Ip));
dvdt = @(I_ext,V,n,m,h,open_k,open_na) (I_ext-Ik(V,n,m,h,open_k)...
    -Ina(V,n,m,h,open_na)-(Il(V,n,m,h)-Ip))/capacitance;
dndt = @(V,n,m,h) (alpha_n(V,n,m,h) * (1-n)) - (beta_n(V,n,m,h)*n);
dmdt = @(V,n,m,h) (alpha_m(V,n,m,h) * (1-m)) - (beta_n(V,n,m,h)*m);
dhdt = @(V,n,m,h) (alpha_h(V,n,m,h) * (1-h)) - (beta_n(V,n,m,h)*h);

% running the simulation from when the external current is applied to the
% end of the simulation 
for i = null_steps:total_steps+1
    % calculating the concentrations inside and outside of the cell for 
    % sodium (Na) and potassium (K)
    Nai = Nai + Ina(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1), ...
        open_val_na) + Il(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1));
    Nao = Nao - Ina(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1),...
        open_val_na) - Il(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1));
    Ki = Ki + Ik(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1), ...
        open_val_k) + Il(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1));
    Ko = Ko - Ik(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1), ...
        open_val_k) - Il(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1));
    % when the external current is stopped the value of the external 
    % current is 0
    if i >= step_off 
        I_ext = 0;
    end
    % when the voltage reaches -55 mV, the sodium gate opens, symbolized by
    % the open_val_na switching from 0 to 1
    if V_vals(i-1) >= na_open && open_val_k == 0
        open_val_na = 1;
    end
    % when the voltage reaches 49.3 mV, the sodium gate closes and 
    % the potassium gate opens, symbolized by the open_val_na switching 
    % from 1 to 0 and the open_val_k switching from 0 to 1
    if V_vals(i-1) >= na_close
        open_val_na = 0;
        open_val_k = 1;
    end
    % when the voltage reaches -77 mV, the potassium gate closes, 
    % symbolized by the open_val_na switching from 1 to 0
    if V_vals(i-1) <= k_close
        open_val_k = 0;
        disp(V_vals(i-1))
    end
    
    % if the sodium concentration inside the cell and potassium
    % concentration outside the cell are positive then the pump is
    % activated, if not the pump is off 
    if Nai > 0 && Ko > 0
        Nai = Nai - Na_flow; % concentration at the previous time step 
        % minus the 3 sodium ions that move out of the cell
        Nao = Nao + Na_flow;  % concentration at the previous time step 
        % plus the 3 sodium ions that move out of the cell
        Ki = Ki - K_flow; % concentration at the previous time step 
        % minus the 2 potassium ions that move into the cell
        Ko = Ko + K_flow; % concentration at the previous time step 
        % plus the 2 potassium ions that into of the cell
        Ip = Il(init_V,init_n,init_m,init_h); % setting the pump constant
    else
        Ip = 0; % if the sodium concentration inside the cell and the 
        % potassium concentration outside the cell are negative then the
        % pump is turned off and the pump constant is set to 0
    end
    
    % Delta 1 uses euler's method to calculate the change in action 
    % potential and gating variables at the given time step
    deltaV_1 = dvdt(I_ext,V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1), open_val_k, open_val_na)* ...
        time_step; 
    deltan_1 = dndt(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1))* ...
    time_step; 
    deltam_1 = dmdt(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1))* ...
        time_step;
    deltah_1 = dhdt(V_vals(i-1),n_vals(i-1),m_vals(i-1),h_vals(i-1))* ...
        time_step;
    % Delta 2 estimates halfway between the previous time step and the
    % current time step using delta 1
    deltaV_2 = dvdt(I_ext,V_vals(i-1)+(deltaV_1*0.5),n_vals(i-1)+(deltan_1*0.5), ...
        m_vals(i-1)+(deltam_1*0.5),h_vals(i-1)+(deltah_1*0.5), open_val_k, ...
        open_val_na)* time_step; 
    deltan_2 = dndt(V_vals(i-1)+(deltaV_1*0.5),n_vals(i-1)+(deltan_1*0.5), ...
        m_vals(i-1)+(deltam_1*0.5),h_vals(i-1)+(deltah_1*0.5))* time_step; 
    deltam_2 = dmdt(V_vals(i-1)+(deltaV_1*0.5),n_vals(i-1)+(deltan_1*0.5), ...
        m_vals(i-1)+(deltam_1*0.5),h_vals(i-1)+(deltah_1*0.5))* time_step;
    deltah_2 = dhdt(V_vals(i-1)+(deltaV_1*0.5),n_vals(i-1)+(deltan_1*0.5),...
        m_vals(i-1)+(deltam_1*0.5),h_vals(i-1)+(deltah_1*0.5))* time_step;
    % Delta 3 estimates halfway between the previous time step and the
    % current time step using delta 2
    deltaV_3 = dvdt(I_ext,V_vals(i-1)+(deltaV_2*0.5),n_vals(i-1)+(deltan_2*0.5), ...
        m_vals(i-1)+(deltam_2*0.5),h_vals(i-1)+(deltah_2*0.5),open_val_k, ...
        open_val_na)* time_step; 
    deltan_3 = dndt(V_vals(i-1)+(deltaV_2*0.5),n_vals(i-1)+(deltan_2*0.5), ...
        m_vals(i-1)+(deltam_2*0.5),h_vals(i-1)+(deltah_2*0.5))* time_step; 
    deltam_3 = dmdt(V_vals(i-1)+(deltaV_2*0.5),n_vals(i-1)+(deltan_2*0.5), ...
        m_vals(i-1)+(deltam_2*0.5),h_vals(i-1)+(deltah_2*0.5))* time_step;
    deltah_3 = dhdt(V_vals(i-1)+(deltaV_2*0.5),n_vals(i-1)+(deltan_2*0.5),...
        m_vals(i-1)+(deltam_2*0.5),h_vals(i-1)+(deltah_2*0.5))* time_step;
    % Delta 4 estimates the change in action potential and gating variables 
    % using delta 3 
    deltaV_4 = dvdt(I_ext,V_vals(i-1)+deltaV_3,n_vals(i-1)+deltan_3, ...
        m_vals(i-1)+deltam_3,h_vals(i-1)+deltah_3,open_val_k, open_val_na)* ...
        time_step; 
    deltan_4 = dndt(V_vals(i-1)+deltaV_3,n_vals(i-1)+deltan_3, ...
        m_vals(i-1)+deltam_3,h_vals(i-1)+deltah_3)* time_step; 
    deltam_4 = dmdt(V_vals(i-1)+deltaV_3,n_vals(i-1)+deltan_3, ...
        m_vals(i-1)+deltam_3,h_vals(i-1)+deltah_3)* time_step;
    deltah_4 = dhdt(V_vals(i-1)+deltaV_3,n_vals(i-1)+deltan_3,...
        m_vals(i-1)+deltam_3,h_vals(i-1)+deltah_3)* time_step;
    % a weighted average is calculated that places more weight on delta 2
    % and delta 3 to calculate the estimated change in action potential and
    % gating variables
    delta_V = (deltaV_1 + (2*deltaV_2) + (2*deltaV_3) + deltaV_4)/6; 
    delta_n = (deltan_1 + (2*deltan_2) + (2*deltan_3) + deltan_4)/6;
    delta_m = (deltam_1 + (2*deltam_2) + (2*deltam_3) + deltam_4)/6;
    delta_h = (deltah_1 + (2*deltah_2) + (2*deltah_3) + deltah_4)/6;
    % the estimated change is added to the previous value of action 
    % potential and gating variables 
    V_vals(i) = V_vals(i-1) + delta_V;
    n_vals(i) = n_vals(i-1) + delta_n;
    m_vals(i) = m_vals(i-1) + delta_m;
    h_vals(i) = h_vals(i-1) + delta_h;
end

% graphing the gating variables with respect to time
figure
plot(0:time_step:simulation_length,n_vals, 0:time_step:simulation_length,...
    m_vals,0:time_step:simulation_length,h_vals)
xlabel("time (ms)")
legend("n","m","h")
title("Gating constants of Hodgkin-Huxley Model with Na-K Pump and Tracking Voltage")

% graphing the action potential with respect to time
figure
plot(0:time_step:simulation_length, V_vals)
xlabel("time (ms)")
title("Action Potential of Hodgkin-Huxley Model with Na-K Pump and Tracking Voltage")
