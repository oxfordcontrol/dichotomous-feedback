clear all
params.volume = 1;%602;

params.kt = 6.12*1e+3;%4.15;
params.kp = 0.00294*60;%1.37; 

params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;


d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 0.167;
params.Rtot   = 6%5.833; 
params.tlat_Taz    = params.delta   * params.Ttot;

params.tx_gfp = 1%.5;

% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;
ops = odeset('AbsTol', 1e-14, 'RelTol', 1e-4);
asp_grid = [logspace(-1, log10(2), 4), linspace(3, 10, 6)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);

x0l = zeros(6*n_t,1);
params.kap_taz = .02.*(params.Atot).'./(1e+3+(params.Atot).');

input_ompr = 6%5.83;
input_p = [0, 0.1, 0.4, 0.7, 1];
T  = 1e+5;
omprp_tc_cl = zeros(length(input_ompr), length(input_p), n_t);
x_tc_cl = zeros(length(input_ompr), length(input_p), n_t);

x_tc_cl_all = cell(length(input_ompr), length(input_p));
t_tc_cl_all = cell(length(input_ompr), length(input_p));

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.P         = input_p(jj);        
        [t_tc_cl_all{ii,jj}, x_tc_cl_all{ii,jj}] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0l, ops);
        omprp_tc_cl(ii,jj,:)  = x_tc_cl_all{ii,jj}(end, [3*n_t + 1 : 4*n_t]);
        x_tc_cl(ii,jj,:)      = x_tc_cl_all{ii,jj}(end, [5*n_t + 1 : 6*n_t]);

    end
end 


input_tazc = [0, logspace(-2, log10(5), 30), linspace(5, 50, 45)];
omprp_tc_ol = zeros(length(input_ompr), length(input_p), n_t);
x_tc_ol = zeros(length(input_ompr), length(input_p), n_t);

x_tc_ol_all = cell(length(input_ompr), length(input_p));
t_tc_ol_all = cell(length(input_ompr), length(input_p));

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_tazc)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.tlat_tazc = params.delta*input_tazc(jj);
        [t_tc_ol_all{ii,jj}, x_tc_ol_all{ii,jj}] = ode15s(@(t,x)two_comp_tazc_ol(t,x,params),[0,T],x0l, ops);
        omprp_tc_ol(ii,jj,:)  = x_tc_ol_all{ii,jj}(end, [3*n_t + 1 : 4*n_t]);
        x_tc_ol(ii,jj,:)      = x_tc_ol_all{ii,jj}(end, [5*n_t + 1 : 6*n_t]);

    end
end

omprp_tc_app = zeros(length(input_ompr), length(input_p), n_t);

for ii = 1:length(input_ompr)
    for jj = 1:length(input_tazc)
        omprp_tc_app(ii,jj,:) = params.kap_taz * params.Ttot/(params.kp*params.Ttot + params.kpc*input_tazc(jj) + params.delta)*params.kt*input_ompr(ii)/(params.kt*input_ompr(ii) + params.delta);
    end
end

x = linspace(0, .1, 100);
num   = (x./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)))/params.delta;

figure(1)
clf; 
hold on
colors_ = {'k', 'c', 'g', 'm', 'r'}

for jj = 2:length(input_p)
    plot(input_p(jj)*tl_x, x, colors_{jj}, 'Linewidth',2);
end
for ii = 1:length(input_ompr)
    for kk = 1:length(params.Atot)
       plot(input_tazc, (squeeze(omprp_tc_ol(ii, :, kk))),'g', 'Linewidth', 3)
       plot(input_tazc, (squeeze(omprp_tc_app(ii, :, kk))),'k--', 'Linewidth', 3)
    end
end
err = (omprp_tc_ol-omprp_tc_app)./omprp_tc_ol;
sqrt(sum(abs(err(:)).^2))
max(abs(err(:)))
ylabel('RR_p Expression')
xlabel('Total PH')
