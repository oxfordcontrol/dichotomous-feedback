clear all
params.volume = 602.2;

params.kt = 6.12*1e+3/params.volume;
params.kp = 4*1e-3/params.volume;
params.ktc = params.kt;
params.kpc = params.kp;

params.Kdr  = 3.14*1e-2*params.volume;
params.hill_coeff = 4;

d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 100;%0.167*params.volume;
params.Rtot   = 3500;%6*params.volume; 
params.tlat_Taz    = params.delta   * params.Ttot;

params.tx_gfp = 1*params.volume;%.5;


ops = odeset('AbsTol', 1e-20, 'RelTol', 1e-7);
% asp_grid = [0, logspace(-2, log10(2),10), linspace(2.01, 6, 10), linspace(6.1, 10, 5)];
asp_grid = [0, logspace(-2, log10(2),40), linspace(2.01, 6, 20), linspace(6.1, 10, 20)];
params.Atot = asp_grid;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);


x0_rr = zeros(7*n_t,1);
x0_tc = zeros(6*n_t,1);

k_a_max = 0.02;%0.1;
K_d_a = 2;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr =1;
n_trc = 20;
input_ompr = params.Rtot;%linspace(.1, 10, n_tr);
input_p    = [0, 0.1, 0.4, 0.7, 1];%linspace(0, 1, n_trc); 0, 0.1, 0.25, 0.5, 0.75, 1
T  = 1e+4;


taz_rr_ss    = zeros(length(input_ompr), length(input_p), n_t);
tazp_rr_ss   = zeros(length(input_ompr), length(input_p), n_t);
ompr_rr_ss   = zeros(length(input_ompr), length(input_p), n_t);
omprp_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprc_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprcp_rr_ss = zeros(length(input_ompr), length(input_p), n_t);
x_prot_rr_ss = zeros(length(input_ompr), length(input_p), n_t);
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.P = input_p(jj);        
        [t_hill, x_hill] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0_rr, ops);
        taz_rr_ss(ii,jj,:)      = x_hill(end,[0*n_t + 1 : 1*n_t]);
        tazp_rr_ss(ii,jj,:)     = x_hill(end,[1*n_t + 1 : 2*n_t]);
        ompr_rr_ss(ii,jj,:)     = x_hill(end,[2*n_t + 1 : 3*n_t]);
        omprp_rr_ss(ii,jj,:)    = x_hill(end,[3*n_t + 1 : 4*n_t]);
        omprc_rr_ss(ii,jj,:)    = x_hill(end,[4*n_t + 1 : 5*n_t]);
        omprcp_rr_ss(ii,jj,:)   = x_hill(end,[5*n_t + 1 : 6*n_t]);
        x_prot_rr_ss(ii,jj,:)   = x_hill(end,[6*n_t + 1 : 7*n_t]);
    end
end


syms Taz TazP TazC OmpRP OmpR OmpRC OmpRCP GFP kap P tlat_ompr beta_r Asp

num_rr   = (OmpRP/params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x_rr  = (params.tx_gfp*(num_rr./(num_rr + 1)));
tl_rc_rr  = P*tl_x_rr;
vx_rr = [params.tlat_Taz; ... 1 
      params.delta*Taz; ...2 
      params.delta*TazP; ...3 
      kap*Taz; .... 4
      params.kt*TazP*OmpR; ... 5
      params.ktc*TazP*OmpRC;... 6
      beta_r;... 7
      params.delta*OmpR;... 8
      params.delta*OmpRP;... 9
      params.kp.*Taz.*OmpRP;... 10
      tl_rc_rr;... 11
      params.delta*OmpRC;... 12
      params.delta*OmpRCP;... 13
      params.kpc.*Taz.*OmpRCP;... 14
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

std_rr    = zeros(length(input_ompr), length(input_p), n_t);
mean_rr   = x_prot_rr_ss;

for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        for kk = 1 : length(params.Atot)
            A_val_rr = eval(subs(Ap_rr, [Taz; TazP; OmpR; OmpRP; OmpRC; OmpRCP; GFP; kap; beta_r; P], ...
                                   [taz_rr_ss(ii,jj,kk);   tazp_rr_ss(ii,jj,kk);   ...
                                    ompr_rr_ss(ii,jj,kk);  omprp_rr_ss(ii,jj,kk); ...
                                    omprc_rr_ss(ii,jj,kk); omprcp_rr_ss(ii,jj,kk); ...
                                    x_prot_rr_ss(ii,jj,kk); ...
                                    params.kap_taz(kk); input_ompr(ii)*params.delta; input_p(jj)]));
            BBt_val_rr = eval(subs(BBt_rr, [Taz; TazP; OmpR; OmpRP; OmpRC; OmpRCP; GFP; kap; beta_r; P], ...
                                     [taz_rr_ss(ii,jj,kk);   tazp_rr_ss(ii,jj,kk);   ...
                                      ompr_rr_ss(ii,jj,kk);  omprp_rr_ss(ii,jj,kk); ...
                                      omprc_rr_ss(ii,jj,kk); omprcp_rr_ss(ii,jj,kk); ...
                                      x_prot_rr_ss(ii,jj,kk); ...
                                      params.kap_taz(kk); input_ompr(ii)*params.delta; input_p(jj)]));                                  
           P_rr = lyap(A_val_rr, BBt_val_rr);
           std_rr(ii,jj,kk)  = sqrt(trace(Cp_rr*P_rr*Cp_rr'));
        end
    end
end
colors_ = {'k', 'c', 'b', 'g', 'r', 'm'}

figure(1)
clf
hold on
ii_x = 1;
plot(squeeze(mean_rr(ii_x,1,:)), squeeze(std_rr(ii_x,1,:)./mean_rr(ii_x,1,:)), 'k', 'Linewidth', 3)
for ii = ii_x%1 : length(input_ompr)
    for jj = 2 : length(input_p)
          plot(squeeze(mean_rr(ii,jj,:)), squeeze(std_rr(ii,jj,:)./mean_rr(ii,jj,:)),strcat(colors_{jj},'-'), 'Linewidth', 3)
    end
end
% axis([0.02, 50, 0.25, 50 ])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Output')
ylabel('Coefficient of Variation')
