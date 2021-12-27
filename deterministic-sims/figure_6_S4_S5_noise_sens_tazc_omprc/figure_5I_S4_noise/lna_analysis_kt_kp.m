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


ops = odeset('AbsTol', 1e-20, 'RelTol', 1e-7);
asp_grid = [0, logspace(-2, log10(2),40), linspace(2.01, 6, 20), linspace(6.1, 30, 20)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);


x0_rr = zeros(7*n_t,1);
x0_tc = zeros(6*n_t,1);

k_a_max = 0.02;%0.1;
K_d_a = 2e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 20;
input_ompr = params.Rtot;%linspace(.1, 10, n_tr);
input_p    = [0,  1];%linspace(0, 1, n_trc);
T  = 1e+4;

params.tlat_ompr = params.delta*params.Rtot;

taz_rr_ss    = zeros(length(kt_grid), length(kp_grid), 2, n_t);
tazp_rr_ss   = zeros(length(kt_grid), length(kp_grid), 2,n_t);
ompr_rr_ss   = zeros(length(kt_grid), length(kp_grid), 2, n_t);
omprp_rr_ss  = zeros(length(kt_grid), length(kp_grid), 2, n_t);
omprc_rr_ss  = zeros(length(kt_grid), length(kp_grid), 2, n_t);
omprcp_rr_ss = zeros(length(kt_grid), length(kp_grid), 2, n_t);
x_prot_rr_ss = zeros(length(kt_grid), length(kp_grid), 2, n_t);
for ii = 1 : length(kt_grid)
    for jj = 1 : length(kp_grid)
        for kk = 1:2
        params.kt  = kt_grid(ii);
        params.kp  = kp_grid(jj);
        params.ktc = params.kt;
        params.kpc = params.kp;
        params.P = input_p(kk);        
        [t_hill, x_hill] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0_rr, ops);
        taz_rr_ss(ii,jj,kk,:)      = x_hill(end,[0*n_t + 1 : 1*n_t]);
        tazp_rr_ss(ii,jj,kk,:)     = x_hill(end,[1*n_t + 1 : 2*n_t]);
        ompr_rr_ss(ii,jj,kk,:)     = x_hill(end,[2*n_t + 1 : 3*n_t]);
        omprp_rr_ss(ii,jj,kk,:)    = x_hill(end,[3*n_t + 1 : 4*n_t]);
        omprc_rr_ss(ii,jj,kk,:)    = x_hill(end,[4*n_t + 1 : 5*n_t]);
        omprcp_rr_ss(ii,jj,kk,:)   = x_hill(end,[5*n_t + 1 : 6*n_t]);
        x_prot_rr_ss(ii,jj,kk,:)   = x_hill(end,[6*n_t + 1 : 7*n_t]);
        end
    end
end


syms Taz TazP TazC OmpRP OmpR OmpRC OmpRCP GFP kap P tlat_ompr beta_r Asp kt ktc kp kpc

num_rr   = (OmpRP/params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x_rr  = (params.tx_gfp*(num_rr./(num_rr + 1)));
tl_rc_rr  = P*tl_x_rr;
vx_rr = [params.tlat_Taz; ... 1 
      params.delta*Taz; ...2 
      params.delta*TazP; ...3 
      kap*Taz; .... 4
      kt*TazP*OmpR; ... 5
      ktc*TazP*OmpRC;... 6
      beta_r;... 7
      params.delta*OmpR;... 8
      params.delta*OmpRP;... 9
      kp.*Taz.*OmpRP;... 10
      tl_rc_rr;... 11
      params.delta*OmpRC;... 12
      params.delta*OmpRCP;... 13
      kpc.*Taz.*OmpRCP;... 14
      tl_x_rr;... 15
      params.delta*GFP]; %16
       %1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
S_rr = [1,-1, 0,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0,-1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0,-1, 0, 1,-1, 0, 1, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 1, 0, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1,-1, 0, 1, 0, 0;...
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0;...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1];  

Ap_rr  = jacobian(simplify(S_rr*vx_rr), [Taz; TazP; OmpR; OmpRP; OmpRC; OmpRCP; GFP]);
BBt_rr = S_rr*diag(vx_rr)*S_rr';
Cp_rr = zeros(1,7); Cp_rr(7) = 1;

std_rr    = zeros(length(kt_grid), length(kp_grid), length(input_p), n_t);
mean_rr   = x_prot_rr_ss;

for ii = 1 : length(kt_grid)
    for jj = 1 : length(kp_grid)
        for kk = 1 : length(params.Atot)
            for ll = 1:2
            A_val_rr = eval(subs(Ap_rr, [Taz; TazP; OmpR; OmpRP; OmpRC; OmpRCP; GFP; kap; beta_r; P;kt; ktc; kp; kpc], ...
                                   [taz_rr_ss(ii,jj,ll, kk);   tazp_rr_ss(ii,jj,ll, kk);   ...
                                    ompr_rr_ss(ii,jj,ll, kk);  omprp_rr_ss(ii,jj,ll, kk); ...
                                    omprc_rr_ss(ii,jj,ll, kk); omprcp_rr_ss(ii,jj,ll, kk); ...
                                    x_prot_rr_ss(ii,jj,ll, kk); ...
                                    params.kap_taz(kk); params.Rtot*params.delta; input_p(ll); ...
                                    kt_grid(ii);kt_grid(ii);kp_grid(jj); kp_grid(jj)]));
            BBt_val_rr = eval(subs(BBt_rr, [Taz; TazP; OmpR; OmpRP; OmpRC; OmpRCP; GFP; kap; beta_r; P;kt; ktc; kp; kpc], ...
                                     [taz_rr_ss(ii,jj,ll, kk);   tazp_rr_ss(ii,jj,ll, kk);   ...
                                      ompr_rr_ss(ii,jj,ll, kk);  omprp_rr_ss(ii,jj,ll, kk); ...
                                      omprc_rr_ss(ii,jj,ll, kk); omprcp_rr_ss(ii,jj,ll, kk); ...
                                      x_prot_rr_ss(ii,jj, ll, kk); ...
                                      params.kap_taz(kk); params.Rtot*params.delta; input_p(ll); ...
                                      kt_grid(ii);kt_grid(ii);kp_grid(jj); kp_grid(jj)]));                                  
           P_rr = lyap(A_val_rr, BBt_val_rr);
           std_rr(ii,jj,ll, kk)  = sqrt(trace(Cp_rr*P_rr*Cp_rr'));
            end
        end
    end
end
colors_ = {'k', 'c', 'b', 'g', 'r'}
style_ = {'',':','--','-.', ':'};

figure(1)
clf
hold on
for ii = 1 : 1%length(kt_grid)
    for jj = 1 : length(kp_grid)
          plot(squeeze(mean_rr(ii, jj, 1,:)), squeeze(std_rr(ii, jj, 1,:)./mean_rr(ii, jj, 1,:)),strcat(colors_{jj},'-'), 'Linewidth', 1)
          plot(squeeze(mean_rr(ii,jj,2, :)), squeeze(std_rr(ii,jj, 2, :)./mean_rr(ii,jj,2, :)),strcat(colors_{jj},'--'))
    end
end
axis([0.02, 50, 0.25, 50 ])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Output')
ylabel('Coefficient of Variation')
