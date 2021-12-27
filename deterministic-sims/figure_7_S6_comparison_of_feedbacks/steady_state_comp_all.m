% clear all
params.volume = 1;%602;
%
kt_grid = [6.12*1e+3, 6.12*1e+1, 1]; %style_  = {'-', '--', ':', '.-'}
kp_grid = [0.00294*60, 4*1e-3, 1*1e-4]; %colors_ = {'k', 'c', 'b', 'g', 'r'}

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;

d_protein =  0.0234;
params.delta = d_protein;
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

k_a_max = 0.02;%0.1;
K_d_a = 2e+3%1e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');
params.kup_hill = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

input_ompr = 6;

% params.kpc = params.kp;
%
% input_omprc = [0, 2, 5, 6, 7];
% xout_tc_ol = zeros(length(input_ompr), length(input_omprc), n_t);
% x_tc_ol = cell(length(input_ompr), length(input_p));
% t_tc_ol = cell(length(input_ompr), length(input_p));
% 
%for ii = 1 : length(input_ompr)
%    for jj = 1 : length(input_omprc)
%        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
%        params.tlat_tazc = params.delta_OmpR*input_omprc(jj);        
%        [t_tc_ol{ii,jj}, x_tc_ol{ii,jj}] = ode15s(@(t,x)two_comp_tazc_ol(t,x,params),[0,T],x0_tc, ops);
%        xout_tc_ol(ii,jj,:) = x_tc_ol{ii,jj}(end,[5*n_t + 1 : 6*n_t]);
%    end
%end
x0_tc = zeros(6*n_t,1);
T = 1e+4;
input_p = [0, 0.1, 0.4, 0.7, 1];
params.tlat_omprc = params.tx_gfp;
xout_tc_cl = zeros(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid), n_t);
x_tc_cl = cell(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid));
t_tc_cl = cell(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid));

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        for kk = 1:length(kt_grid)
            for ll = 1:length(kp_grid)            
                params.kt  = kt_grid(kk);%.04 %6.12*1e-3*1e+3/2;%4.15;
                params.kp  = kp_grid(ll); %1e-3; %0.00294*60*1e-6*1e+3*1e+1*1;%1.37; 
                params.ktc = params.kt;
                params.kpc = params.kp;
                params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
                params.P = input_p(jj);
                [t_tc_cl{ii,jj,kk,ll}, x_tc_cl{ii,jj,kk,ll}] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0_tc, ops);
                xout_tc_cl(ii,jj,kk,ll,:) = x_tc_cl{ii,jj,kk,ll}(end,[5*n_t + 1 : 6*n_t]);
            end
        end
    end
end


% input_omprc = [0, 2, 5, 6, 7];
% xout_rr_ol = zeros(length(input_ompr), length(input_omprc), n_t);
% x_rr_ol = cell(length(input_ompr), length(input_p));
% t_rr_ol = cell(length(input_ompr), length(input_p));
%for ii = 1 : length(input_ompr)
%    for jj = 1 : length(input_omprc)
%        params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
%        params.tlat_omprc = params.delta_OmpR*input_omprc(jj);        
%        [t_rr_ol{ii,jj}, x_rr_ol{ii,jj}] = ode15s(@(t,x)two_comp_sink_hill_ol(t,x,params),[0,T],x0_rr, ops);
%        xout_rr_ol(ii,jj,:) = x_rr_ol{ii,jj}(end,[6*n_t + 1 : 7*n_t]);
%    end
%end

x0_rr = zeros(7*n_t,1);
params.tlat_omprc = params.tx_gfp;
xout_rr_cl = zeros(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid), n_t);
x_rr_cl = cell(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid));
t_rr_cl = cell(length(input_ompr), length(input_p), length(kt_grid), length(kp_grid));

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        for kk = 1:length(kt_grid)
            for ll = 1:length(kp_grid)            
                params.tlat_ompr = params.delta_OmpR*input_ompr(ii);
                params.P   = input_p(jj);
                params.kt  = kt_grid(kk);
                params.kp  = kp_grid(ll);
                params.ktc = params.kt;
                params.kpc = params.kp;                
                [t_rr_cl{ii,jj,kk,ll}, x_rr_cl{ii,jj,kk,ll}] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0_rr, ops);
                xout_rr_cl(ii,jj,kk,ll,:) = x_rr_cl{ii,jj,kk,ll}(end,[6*n_t + 1 : 7*n_t]);
            end
        end
    end
end


max_gfp = max(xout_rr_cl(:));
colors_ = {'k', 'c', 'b', 'g', 'r'}
style_  = {'-', '--', ':', '.-'}

figure(1)
clf; 
hold on
kk = 2;
ll = 2;
for ii = 1:length(input_ompr)
    for jj = 1:length(input_p)
       plot(asp_grid, squeeze(xout_tc_cl(ii, jj, kk, ll,:))/max_gfp, strcat(colors_{jj}, '--'), 'LineWidth',2)
       plot(asp_grid, squeeze(xout_rr_cl(ii, jj, kk, ll,:))/max_gfp, strcat(colors_{jj}, '-'), 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')


figure(2)
clf; 
hold on
ii = 1;
jj = 5;
for kk = 1: 1%length(kt_grid)
    for ll = 1:length(kp_grid)
       col_tc = strcat(colors_{ll}, style_{kk});
       col_rr = strcat(colors_{ll}, style_{kk+1});
       plot(asp_grid, squeeze(xout_tc_cl(ii, 1, kk,ll,:))/max_gfp, strcat(colors_{ll}, '-'), 'LineWidth',1)
     %  plot(asp_grid, squeeze(xout_tc_cl(ii,jj,kk,ll,:))/max_gfp, col_tc, 'LineWidth',1)
       plot(asp_grid, squeeze(xout_rr_cl(ii,jj,kk,ll,:))/max_gfp, col_rr, 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')
figure(3)
clf; 
hold on
ii = 1;
jj = 5;
for kk = 1: 1%length(kt_grid)
    for ll = 1:length(kp_grid)
       col_tc = strcat(colors_{ll}, style_{kk+1});
       col_rr = strcat(colors_{ll}, strcat('x',style_{kk+1}));
       plot(asp_grid, squeeze(xout_tc_cl(ii, 1, kk,ll,:))/max_gfp, strcat(colors_{ll}, '-'), 'LineWidth',1)
       plot(asp_grid, squeeze(xout_tc_cl(ii,jj,kk,ll,:))/max_gfp, col_tc, 'LineWidth',2)
      % plot(asp_grid, squeeze(xout_rr_cl(ii,jj,kk,ll,:))/max_gfp, col_rr, 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')
figure(4)
clf; 
hold on
ii = 1;
jj = 2;
for kk = 1:length(kt_grid)
    for ll = 1:length(kp_grid)
       col_tc = strcat(colors_{ll}, style_{kk});
       col_rr = strcat(colors_{ll}, strcat('x',style_{kk}));
%        plot(asp_grid, squeeze(xout_tc_cl(ii, 1, kk,ll,:))/max_gfp, 'k-', 'LineWidth',4)
       plot(asp_grid, squeeze(xout_tc_cl(ii,jj,kk,ll,:))/max_gfp, col_tc, 'LineWidth',1)
%        plot(asp_grid, squeeze(xout_rr_cl(ii,jj,kk,ll,:))/max_gfp, col_rr, 'LineWidth',2)
    end
end
ylabel('Output')
xlabel('Inducer')
