% clear all
params.volume = 1;%602;

params.kt = 6.12*1e+3;%4.15;
params.kp = 4*1e-2;%1.37; 
params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;

d_protein =  0.0234;
params.delta_K    = d_protein; %0.03;
params.delta_OmpR  = d_protein; %0.03;
params.delta_GFP  = d_protein; %0.03;
params.delta_Taz = d_protein;
params.delta_Asp  = d_protein;


params.Ttot   = 0.167;
params.tlat_Taz    = params.delta_Taz   * params.Ttot;
% params.tlat_OmpR   = params.delta_OmpR  * params.Rtot;

params.tx_gfp = 1;%.5;


% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;
ops = odeset('AbsTol', 1e-14, 'RelTol', 1e-4);
asp_grid = [0, logspace(-1, log10(2),20), linspace(2.01, 10, 20)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);

x0l = zeros(6*n_t,1);
k_a_max = 0.02;%0.1;
K_d_a = 2e+3;
params.kup_hill = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 5;
input_ompr = 6;

x0_rr = zeros(8*n_t,1);
params.kpc = params.kp;

T= 1e+4

input_p = [0, 0.1, 0.4, 0.7, 1];
input_q = [0, 0.1, 0.4, 0.7, 1];
x_prot_cl = zeros(length(input_ompr), length(input_p), n_t);
x_rr_cl = cell(length(input_ompr), length(input_p));
t_rr_cl = cell(length(input_ompr), length(input_p));


params.tlat_omprc = params.tx_gfp;
for ii = 1 : length(input_q)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta_OmpR*input_ompr;
        params.P = input_p(jj);
        params.Q = input_p(ii);
        [t_rr_cl{ii,jj}, x_rr_cl{ii,jj}] = ode15s(@(t,x)two_comp_sr_ph_feed(t,x,params),[0,T],x0_rr, ops);
        x_prot_cl(ii,jj,:) = x_rr_cl{ii,jj}(end,[7*n_t + 1 : 8*n_t]);
    end
end


max_gfp = max(x_prot_cl(:));
colors_ = {'k', 'c', 'b', 'g', 'r'}
style_ = {'-', ':', '--', '.-', 'x-'}

figure(1)
clf; 
hold on
for ii = 1:length(input_q)
    for jj = 1:length(input_p)
       plot(asp_grid, squeeze(x_prot_cl(ii,jj,:))/max_gfp, strcat(colors_{jj}, style_{ii}), 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')

