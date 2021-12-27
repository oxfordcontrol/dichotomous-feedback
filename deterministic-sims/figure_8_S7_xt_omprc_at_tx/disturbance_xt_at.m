clear all
params.volume = 1;%602;

params.kup_xt = 0.005%6*1e-3;
params.kt = 6.12*1e+3%4.15;
params.kp = 4*1e-3;%1.37; 


params.kpc = 6.12*1e+3/1.5*1e-6;% 0.004;%params.kp;
params.ktx     = params.kt;
params.kpx     = params.kp;
params.kp_tazc = params.kp;

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;

d_protein = 0.0234;
params.delta    = d_protein; %0.03;
params.delta_OmpR  = d_protein; %0.03;
params.delta_GFP  = d_protein; %0.03;
params.delta_Taz = d_protein;
params.delta_Asp  = d_protein;


params.Ttot   = 0.167;
params.Rtot   = 6; 
params.tlat_Taz    = params.delta_Taz   * params.Ttot;

params.tx_gfp = 1;%.5;


% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;

input_ompr = params.Rtot;
input_p = [0, 0.1, 0.4, 0.7,  1];

asp_grid = [0, logspace(-2, 1, 20)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);

k_a_max = 0.02;%0.1;
K_d_a = 2e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

T= 2000;
x0_rr = zeros(8*n_t,1);
params.tlat_ompr = params.delta_OmpR*input_ompr(1);

dist_signals = @(t)([0; 0; 0; 0; 0; 0]);

ops = odeset('AbsTol', 1e-14, 'RelTol', 1e-10);
params.ktx     = params.kt;
params.kpx     = params.kp;
params.k_titr = 0;
params.P = 0;
[t_rr_cl_0,  x_rr_cl_0]    = ode15s(@(t,x)two_comp_taz_at(t, x, params), [0,T], x0_rr, ops);
params.P = 0.1;
[t_rr_cl_01,  x_rr_cl_01]  = ode15s(@(t,x)two_comp_taz_at(t, x, params), [0,T], x0_rr, ops);
params.P = 0.4;
[t_rr_cl_04,  x_rr_cl_04]  = ode15s(@(t,x)two_comp_taz_at(t, x, params), [0,T], x0_rr, ops);
params.P = 0.7;
[t_rr_cl_07,  x_rr_cl_07]  = ode15s(@(t,x)two_comp_taz_at(t, x, params), [0,T], x0_rr, ops);
params.P = 1;
[t_rr_cl_1,  x_rr_cl_1]    = ode15s(@(t,x)two_comp_taz_at(t, x, params), [0,T], x0_rr, ops);


colors_ = {'k', 'c', 'b', 'g', 'r'}

figure(1)
clf;
hold on
index_x  = 5*n_t+1:6*n_t;
max_ = max(x_rr_cl_0(end, index_x));
plot(asp_grid, x_rr_cl_0(end, index_x)/max_,colors_{1})
plot(asp_grid, x_rr_cl_01(end, index_x)/max_,colors_{2})
plot(asp_grid, x_rr_cl_04(end, index_x)/max_,colors_{3})
plot(asp_grid, x_rr_cl_07(end, index_x)/max_,colors_{4})
plot(asp_grid, x_rr_cl_1(end, index_x)/max_,colors_{5})

figure(2)
clf;
hold on
index_x  = 5*n_t+1:6*n_t;
plot(asp_grid, x_rr_cl_0(end, index_x)-x_rr_cl_0(end, 5*n_t+1),colors_{1})
plot(asp_grid, x_rr_cl_01(end, index_x)-x_rr_cl_01(end, 5*n_t+1),colors_{2})
plot(asp_grid, x_rr_cl_04(end, index_x)-x_rr_cl_04(end, 5*n_t+1),colors_{3})
plot(asp_grid, x_rr_cl_07(end, index_x)-x_rr_cl_07(end, 5*n_t+1),colors_{4})
plot(asp_grid, x_rr_cl_1(end, index_x)-x_rr_cl_1(end, 5*n_t+1),colors_{5})

figure(3)
clf;
hold on
data_bar = [x_rr_cl_0(end,5*n_t+1), x_rr_cl_01(end,5*n_t+1), x_rr_cl_04(end,5*n_t+1), x_rr_cl_07(end,5*n_t+1), x_rr_cl_1(end,5*n_t+1)]/max_;
barh(data_bar)
