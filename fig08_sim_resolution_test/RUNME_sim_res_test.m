%% Run CToF-DOT: Fig. 8

rng(1);

od_vec = [0]; % Factor determining amount of excitation light reduced by an emission filter
                  
sim_res_test;

% Run CToF-DOT: Supplemental Material

od_vec = [5]; % Factor determining amount of excitation light reduced by an emission filter
                  
sim_res_test;

%% Run DOT Fig. 8

rng(1);

od_vec = [0]; % Factor determining amount of excitation light reduced by an emission filter
                  
sim_dot_res_test;

% Run DOT Supplemental Material

od_vec = [5]; % Factor determining amount of excitation light reduced by an emission filter
                  
sim_dot_res_test;
