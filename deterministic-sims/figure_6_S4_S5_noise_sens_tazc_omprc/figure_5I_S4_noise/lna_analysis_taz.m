clear all
params.volume = 602.2;

params.kt = 6.12*1e+3/params.volume;
params.kp = 4*1e-3/params.volume;
params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2*params.volume;
params.hill_coeff = 2;

d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 100;%0.167;
params.Rtot   = 3500;%5.833; 
params.tlat_Taz    = params.delta   * params.Ttot;

params.tx_gfp = 1*params.volume;%.5;


% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;
ops = odeset('AbsTol', 1e-20, 'RelTol', 1e-7);
asp_grid = [0, logspace(-2, log10(2),20), linspace(2.1, 30, 20)];
params.Atot = asp_grid;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);


x0_rr = zeros(7*n_t,1);
x0_tc = zeros(6*n_t,1);

k_a_max = 0.02;
K_d_a = 2;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 20;
input_ompr = params.Rtot;%linspace(.1, 10, n_tr);
input_p    = [0, 0.1, 0.4, 0.7, 1];%linspace(0, 1, n_trc);
T  = 1e+4;
taz_tc_ss    = zeros(length(input_ompr), length(input_p), n_t);
tazp_tc_ss   = zeros(length(input_ompr), length(input_p), n_t);
ompr_tc_ss   = zeros(length(input_ompr), length(input_p), n_t);
omprp_tc_ss  = zeros(length(input_ompr), length(input_p), n_t);
tazc_tc_ss  = zeros(length(input_ompr), length(input_p), n_t);
x_prot_tc_ss = zeros(length(input_ompr), length(input_p), n_t);
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.P = input_p(jj);        
        [t_hill, x_hill] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0_tc, ops);
        taz_tc_ss(ii,jj,:)      = x_hill(end,[0*n_t + 1 : 1*n_t]);
        tazp_tc_ss(ii,jj,:)     = x_hill(end,[1*n_t + 1 : 2*n_t]);
        ompr_tc_ss(ii,jj,:)     = x_hill(end,[2*n_t + 1 : 3*n_t]);
        omprp_tc_ss(ii,jj,:)    = x_hill(end,[3*n_t + 1 : 4*n_t]);
        tazc_tc_ss(ii,jj,:)     = x_hill(end,[4*n_t + 1 : 5*n_t]);
        x_prot_tc_ss(ii,jj,:)   = x_hill(end,[5*n_t + 1 : 6*n_t]);
    end
end



syms Taz TazP TazC OmpRP OmpR OmpRC OmpRCP GFP kap P tlat_ompr beta_r Asp

num_tc   = (OmpRP/params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x_tc  = (params.tx_gfp*(num_tc./(num_tc + 1)));
tl_tc  = P*tl_x_tc;
vx_tc = [params.tlat_Taz; ... 1 
      params.delta*Taz; ...2 
      params.delta*TazP; ...3 
      kap*Taz; .... 4
      params.kt*TazP*OmpR; ... 5
      beta_r;... 6
      params.delta*OmpR;... 7
      params.delta*OmpRP;... 8
      params.kp.*Taz.*OmpRP;... 9
      tl_tc;... 10
      params.delta*TazC;... 11
      params.kpc.*TazC.*OmpRP;... 12
      tl_x_tc;... 13
      params.delta*GFP]; %14
    %   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
S_tc = [1,-1, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0;...
        0, 0, 0, 0, 1, 0, 0,-1,-1, 0, 0,-1, 0, 0;...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0;...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1];  

  
Ap_tc  = jacobian(simplify(S_tc*vx_tc), [Taz; TazP; OmpR; OmpRP; TazC; GFP]);
BBt_tc = S_tc*diag(vx_tc)*S_tc';
Cp_tc = zeros(1,6); Cp_tc(6) = 1;

std_tc    = zeros(length(input_ompr), length(input_p), n_t);
mean_tc   = x_prot_tc_ss;

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        for kk = 1 : length(params.Atot)
           A_val_tc = eval(subs(Ap_tc, [Taz; TazP; OmpR; OmpRP; TazC; GFP; kap; beta_r; P], ...
                                   [taz_tc_ss(ii,jj,kk);   tazp_tc_ss(ii,jj,kk);   ...
                                    ompr_tc_ss(ii,jj,kk);  omprp_tc_ss(ii,jj,kk); ...
                                    tazc_tc_ss(ii,jj,kk);  x_prot_tc_ss(ii,jj,kk); ...
                                    params.kap_taz(kk); input_ompr(ii)*params.delta; input_p(jj)]));
            BBt_val_tc = eval(subs(BBt_tc, [Taz; TazP; OmpR; OmpRP; TazC;  GFP; kap; beta_r; P], ...
                                     [taz_tc_ss(ii,jj,kk);   tazp_tc_ss(ii,jj,kk);   ...
                                      ompr_tc_ss(ii,jj,kk);  omprp_tc_ss(ii,jj,kk); ...
                                      tazc_tc_ss(ii,jj,kk);  x_prot_tc_ss(ii,jj,kk); ...
                                      params.kap_taz(kk); input_ompr(ii)*params.delta; input_p(jj)]));
           P_tc = lyap(A_val_tc, BBt_val_tc);
           std_tc(ii,jj,kk)  = sqrt(trace(Cp_tc*P_tc*Cp_tc'));
%            mean_(ii,jj,kk) = x_prot_hill_ss(ii,jj,kk);
        end
    end
end
colors_ = {'k', 'c', 'g', 'm', 'r'}

figure(1)
clf
hold on
ii_x = 1;
plot(squeeze(mean_tc(ii_x,1,:)), squeeze(std_tc(ii_x,1,:)./mean_tc(ii_x,1,:)), 'k', 'Linewidth', 3)
for ii = ii_x%1 : length(input_ompr)
    for jj = 2 : length(input_p)
          plot(squeeze(mean_tc(ii,jj,:)), squeeze(std_tc(ii,jj,:)./mean_tc(ii,jj,:)),strcat(colors_{jj},'-.'))
    end
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Output')
ylabel('Coefficient of Variation')
figure(2)
clf
hold on
ii_x = 1;
plot(asp_grid, squeeze(std_tc(ii_x,1,:)./mean_tc(ii_x,1,:)), 'k',  'Linewidth', 3)
for ii = ii_x%1 : length(input_ompr)
    for jj = 2 : length(input_p)
          plot(asp_grid, squeeze(std_tc(ii,jj,:)./mean_tc(ii,jj,:)),strcat(colors_{jj},'-.'))
    end
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Output')
ylabel('Coefficient of Variation')
save lna_data_tazc.mat mean_tc std_tc input_p asp_grid
