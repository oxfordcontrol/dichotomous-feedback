clear all
params.volume = 1;%602;

kt_grid = [6.12*1e+3, 6.12*1e+1, 1]; %style_  = {'-', '--', ':', '.-'}
kp_grid = [0.00294*60, 4*1e-3, 5*1e-2]; %colors_ = {'k', 'c', 'b', 'g', 'r'}

params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;

d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 0.167;
params.Rtot   = 6; 
params.tlat_Taz    = params.delta   * params.Ttot;

params.tx_gfp = 1;%.5;


% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;
ops = odeset('AbsTol', 1e-20, 'RelTol', 1e-7);
asp_grid = [0, logspace(-2, log10(2),20), linspace(2.1, 30, 20)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);


x0_rr = zeros(7*n_t,1);
x0_tc = zeros(6*n_t,1);

k_a_max = 0.02;
K_d_a = 2e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 20;
T  = 1e+4;
params.tlat_ompr = params.delta*params.Rtot;
input_p = [0, 1];
taz_tc_ss    = zeros(length(kt_grid), length(kp_grid), 2, n_t);
tazp_tc_ss   = zeros(length(kt_grid), length(kp_grid), 2, n_t);
ompr_tc_ss   = zeros(length(kt_grid), length(kp_grid), 2, n_t);
omprp_tc_ss  = zeros(length(kt_grid), length(kp_grid), 2, n_t);
tazc_tc_ss   = zeros(length(kt_grid), length(kp_grid), 2, n_t);
x_prot_tc_ss = zeros(length(kt_grid), length(kp_grid), 2, n_t);
for ii = 1 : length(kt_grid)
    for jj = 1 : length(kp_grid)
        for kk = 1:2
        params.kt  = kt_grid(ii);
        params.kp  = kp_grid(jj);

        params.kpc = params.kp;
        params.P = input_p(kk);        
        [t_hill, x_hill] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0_tc, ops);
        taz_tc_ss(ii,jj, kk,:)      = x_hill(end,[0*n_t + 1 : 1*n_t]);
        tazp_tc_ss(ii,jj, kk,:)     = x_hill(end,[1*n_t + 1 : 2*n_t]);
        ompr_tc_ss(ii,jj, kk,:)     = x_hill(end,[2*n_t + 1 : 3*n_t]);
        omprp_tc_ss(ii,jj, kk,:)    = x_hill(end,[3*n_t + 1 : 4*n_t]);
        tazc_tc_ss(ii,jj, kk,:)     = x_hill(end,[4*n_t + 1 : 5*n_t]);
        x_prot_tc_ss(ii,jj, kk,:)   = x_hill(end,[5*n_t + 1 : 6*n_t]);
        end
    end
end



syms Taz TazP TazC OmpRP OmpR OmpRC OmpRCP GFP kap P tlat_ompr beta_r Asp kt kp kpc

num_tc   = (OmpRP/params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x_tc  = (params.tx_gfp*(num_tc./(num_tc + 1)));
tl_tc  = P*tl_x_tc;
vx_tc = [params.tlat_Taz; ... 1 
      params.delta*Taz; ...2 
      params.delta*TazP; ...3 
      kap*Taz; .... 4
      kt*TazP*OmpR; ... 5
      beta_r;... 6
      params.delta*OmpR;... 7
      params.delta*OmpRP;... 8
      kp.*Taz.*OmpRP;... 9
      tl_tc;... 10
      params.delta*TazC;... 11
      kpc.*TazC.*OmpRP;... 12
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

std_tc    = zeros(length(kt_grid), length(kp_grid), n_t);
mean_tc   = x_prot_tc_ss;

for ii = 1 : length(kt_grid)
    for jj = 1 : length(kp_grid)
        for kk = 1 : length(params.Atot)
            for ll = 1:2            
           A_val_tc = eval(subs(Ap_tc, [Taz; TazP; OmpR; OmpRP; TazC; GFP; kap; beta_r; P; kt; kp; kpc], ...
                                   [taz_tc_ss(ii,jj,ll, kk);   tazp_tc_ss(ii,jj,ll, kk);   ...
                                    ompr_tc_ss(ii,jj,ll, kk);  omprp_tc_ss(ii,jj,ll, kk); ...
                                    tazc_tc_ss(ii,jj,ll, kk);  x_prot_tc_ss(ii,jj,ll, kk); ...
                                    params.kap_taz(kk); params.Rtot*params.delta; input_p(ll); ...
                                    kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
            BBt_val_tc = eval(subs(BBt_tc, [Taz; TazP; OmpR; OmpRP; TazC;  GFP; kap; beta_r; P; kt; kp; kpc], ...
                                     [taz_tc_ss(ii,jj, ll, kk);   tazp_tc_ss(ii,jj,ll, kk);   ...
                                      ompr_tc_ss(ii,jj,ll, kk);  omprp_tc_ss(ii,jj,ll, kk); ...
                                      tazc_tc_ss(ii,jj,ll, kk);  x_prot_tc_ss(ii,jj,ll, kk); ...
                                      params.kap_taz(kk); params.Rtot*params.delta; input_p(ll);...
                                      kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
           P_tc = lyap(A_val_tc, BBt_val_tc);
           std_tc(ii,jj, ll, kk)  = sqrt(trace(Cp_tc*P_tc*Cp_tc'));
            end
        end
    end
end
colors_ = {'k', 'c', 'b', 'm', 'r'}
style_ = {'',':','--','-.', ':'};

figure(1)
clf
hold on
for ii = 1 : 1%length(kt_grid)
    for jj = 1 : length(kp_grid)
          plot(squeeze(mean_tc(ii, jj, 1,:)), squeeze(std_tc(ii, jj, 1,:)./mean_tc(ii, jj, 1,:)),strcat(colors_{jj},'-'), 'Linewidth', 1)
          plot(squeeze(mean_tc(ii,jj,2, :)), squeeze(std_tc(ii,jj, 2, :)./mean_tc(ii,jj,2, :)),strcat(colors_{jj},'--') )
    end
end
axis([0.02, 50, 0.25, 50 ])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Output')
ylabel('Coefficient of Variation')
%figure(2)
%clf
%hold on
%ii_x = 1;
%plot(asp_grid, squeeze(std_tc(ii_x,1,1,:)./mean_tc(ii_x,1,1,:)), 'k',  'Linewidth', 3)
%for ii = 1 : length(kt_grid)
%    for jj = 1 : length(kp_grid)
%          plot(asp_grid, squeeze(std_tc(ii,jj,1,:)./mean_tc(ii,jj,1,:)),strcat(colors_{jj},style_{ii}), 'Linewidth', 3)
%          plot(asp_grid, squeeze(std_tc(ii,jj,2,:)./mean_tc(ii,jj,2,:)),strcat(colors_{jj},style_{ii}))
%          
%    end
%end
%
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
