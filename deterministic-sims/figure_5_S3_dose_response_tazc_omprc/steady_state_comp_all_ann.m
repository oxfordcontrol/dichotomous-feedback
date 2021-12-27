% clear all
params.volume = 1;%602;

params.kt = .04 %6.12*1e-3*1e+3/2;%4.15;
params.kp = 1e-3; %0.00294*60*1e-6*1e+3*1e+1*1;%1.37; 
params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 6;

d_protein =  0.0234;
params.delta_K    = d_protein; %0.03;
params.delta   = d_protein; %0.03;
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

k_a_max = 0.02;%0.1;
K_d_a = 2e+3%1e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');
params.kup_hill = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 5;
input_ompr = linspace(1, 10, 10);

x0_tc = zeros(6*n_t,1);
params.kpc = params.kp;
input_p = [0, 0.1, 0.4, 0.7, 1];
T = 1e+4;

xout_tc_cl = zeros(length(input_ompr), length(input_p), n_t);
x_tc_cl = cell(length(input_ompr), length(input_p));
t_tc_cl = cell(length(input_ompr), length(input_p));

params.tlat_omprc = params.tx_gfp;
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
        params.P = input_p(jj);
        [t_tc_cl{ii,jj}, x_tc_cl{ii,jj}] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0_tc, ops);
        xout_tc_cl(ii,jj,:) = x_tc_cl{ii,jj}(end,[5*n_t + 1 : 6*n_t]);
    end
end

x0_tc = zeros(6*n_t,1);
params.kpc = params.kp;

xout_ta_cl = zeros(length(input_ompr), length(input_p), n_t);
x_ta_cl = cell(length(input_ompr), length(input_p));
t_ta_cl = cell(length(input_ompr), length(input_p));

input_p = [0, 0.1, 0.4, 0.7, 1];
params.tlat_omprc = params.tx_gfp;
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
        params.P = input_p(jj);
        [t_ta_cl{ii,jj}, x_ta_cl{ii,jj}] = ode15s(@(t,x)two_comp_taz_at(t,x,params),[0,T],x0_tc, ops);
        xout_ta_cl(ii,jj,:) = x_ta_cl{ii,jj}(end,[5*n_t + 1 : 6*n_t]);
    end
end


x0_tc = zeros(6*n_t,1);
params.kpc = params.kt;

xout_ra_cl = zeros(length(input_ompr), length(input_p), n_t);
x_ra_cl = cell(length(input_ompr), length(input_p));
t_ra_cl = cell(length(input_ompr), length(input_p));

input_p = [0, 0.1, 0.4, 0.7, 1];
params.tlat_omprc = params.tx_gfp;
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
        params.P = input_p(jj);
        [t_ra_cl{ii,jj}, x_ra_cl{ii,jj}] = ode15s(@(t,x)two_comp_ompr_at(t,x,params),[0,T],x0_tc, ops);
        xout_ra_cl(ii,jj,:) = x_ra_cl{ii,jj}(end,[5*n_t + 1 : 6*n_t]);
    end
end


x0_rr = zeros(7*n_t,1);
params.kpc = params.kp;

xout_rr_cl = zeros(length(input_ompr), length(input_p), n_t);
x_rr_cl = cell(length(input_ompr), length(input_p));
t_rr_cl = cell(length(input_ompr), length(input_p));

input_p = [0, 0.1, 0.4, 0.7, 1];
params.tlat_omprc = params.tx_gfp;
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
        params.P = input_p(jj);
        [t_rr_cl{ii,jj}, x_rr_cl{ii,jj}] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0_rr, ops);
        xout_rr_cl(ii,jj,:) = x_rr_cl{ii,jj}(end,[6*n_t + 1 : 7*n_t]);
    end
end


max_gfp = 1;%max(xout_rr_cl(:));
colors_ = {'k', 'c', 'b', 'g', 'r'}

figure(1)
clf; 
hold on
for ii = 5:5%length(input_ompr)
    for jj = 1:length(input_omprc)
       plot(asp_grid, squeeze(xout_ra_cl(ii,jj,:))/max_gfp, strcat(colors_{jj}, ':'), 'LineWidth',4)
       plot(asp_grid, squeeze(xout_tc_cl(ii,jj,:))/max_gfp, strcat(colors_{jj}, '-.'), 'LineWidth',2)
       plot(asp_grid, squeeze(xout_ta_cl(ii,jj,:))/max_gfp, strcat(colors_{jj}, '--'), 'LineWidth',1)
       plot(asp_grid, squeeze(xout_rr_cl(ii,jj,:))/max_gfp, strcat(colors_{jj}, '-'), 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')

figure(2)
clf; 
hold on
for kk = length(asp_grid):length(asp_grid)
    for jj = 1:length(input_omprc)
       plot(input_ompr, squeeze(xout_ra_cl(:,jj,kk))/max_gfp, strcat(colors_{jj}, ':'), 'LineWidth',4)
       plot(input_ompr, squeeze(xout_tc_cl(:,jj,kk))/max_gfp, strcat(colors_{jj}, '-.'), 'LineWidth',2)
       plot(input_ompr, squeeze(xout_ta_cl(:,jj,kk))/max_gfp, strcat(colors_{jj}, '--'), 'LineWidth',1)
       plot(input_ompr, squeeze(xout_rr_cl(:,jj,kk))/max_gfp, strcat(colors_{jj}, '-'), 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')



