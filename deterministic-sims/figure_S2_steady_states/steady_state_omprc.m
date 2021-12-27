clear all
params.volume = 1;%602;

params.kt = 6.12*1e+3;%4.15;
params.kp = 0.00294*60;%1.37; 
%params.kp = 0.004;%1.37; AP changed 16.06.21
params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;


d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 0.167;
params.Rtot   = 6; 
params.tlat_Taz    = params.delta   * params.Ttot;


params.tx_gfp = 1;%.5;


ops = odeset('AbsTol', 1e-14, 'RelTol', 1e-4);
asp_grid = [logspace(-1, log10(2), 4), linspace(3, 10, 6)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);

x0l = zeros(7*n_t,1);
params.kap_taz = .02.*(params.Atot).'./(2e+3+(params.Atot).');

input_ompr = 6;
input_p = [0, 0.1, 0.4, 0.7, 1];
T  = 1e+4;
% omprp_rc_cl = zeros(length(input_ompr), length(input_p), n_t);
% x_rc_cl = zeros(length(input_ompr), length(input_p), n_t);
% 
% x_rc_cl_all = cell(length(input_ompr), length(input_p));
% t_rc_cl_all = cell(length(input_ompr), length(input_p));
% 
% for ii = 1 : length(input_ompr)
%     for jj = 1 : length(input_p)
%         params.tlat_ompr = params.delta*input_ompr(ii);
%         params.P         = input_p(jj);        
%         [t_rc_cl_all{ii,jj}, x_rc_cl_all{ii,jj}] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0l, ops);
%         omprp_rc_cl(ii,jj,:)  = x_rc_cl_all{ii,jj}(end, [3*n_t + 1 : 4*n_t]);
%         x_rc_cl(ii,jj,:)      = x_rc_cl_all{ii,jj}(end, [5*n_t + 1 : 6*n_t]);
% 
%     end
% end 


input_omprc = linspace(0, 20, 40);
omprp_rc_ol = zeros(length(input_ompr), length(input_p), n_t);
x_rc_ol = zeros(length(input_ompr), length(input_p), n_t);

x_rc_ol_all = cell(length(input_ompr), length(input_p));
t_rc_ol_all = cell(length(input_ompr), length(input_p));

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_omprc)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.tlat_omprc = params.delta*input_omprc(jj);
        [t_rc_ol_all{ii,jj}, x_rc_ol_all{ii,jj}] = ode15s(@(t,x)two_comp_sink_hill_ol(t,x,params),[0,T],x0l, ops);
        omprp_rc_ol(ii,jj,:)  = x_rc_ol_all{ii,jj}(end, [3*n_t + 1 : 4*n_t]);
        x_rc_ol(ii,jj,:)      = x_rc_ol_all{ii,jj}(end, [5*n_t + 1 : 6*n_t]);
    end
end

omprp_app = zeros(length(input_ompr), length(input_p), n_t);

for ii = 1:length(input_ompr)
    for jj = 1:length(input_omprc)
        omprp_app(ii,jj,:) = params.kap_taz* params.Ttot/(params.kp*params.Ttot + params.delta) * ...                      
                             params.kt*input_ompr(ii)/(params.kt*input_ompr(ii) + params.ktc*input_omprc(jj)+params.delta);
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
       plot(input_omprc, (squeeze(omprp_rc_ol(ii, :, kk))),'g', 'Linewidth', 3)
       plot(input_omprc, (squeeze(omprp_app(ii, :, kk))),'k--', 'Linewidth', 3)
    end
end
err = (omprp_rc_ol-omprp_app)./omprp_rc_ol;
sqrt(sum(abs(err(:)).^2))
max(abs(err(:)))

ylabel('RR_p Expression')
xlabel('Total SR')
