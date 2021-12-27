
clear all
params.volume = 1;%602;

kt_grid = [6.12*1e+3, 6.12*1e+1, 1];%4.15;
kp_grid = [0.00294*60, 4*1e-2, 4*1e-3];%1.37; 


params.Kdr  = 3.14*1e-2;
params.hill_coeff = 4;

d_protein = 0.0234;
params.delta    = d_protein; %0.03;

params.Ttot   = 0.167;
params.tlat_Taz    = params.delta   * params.Ttot;

% params.tlat_OmpR   = params.delta  * params.Rtot;

params.tx_gfp = 1;%.5;


% params.tx_ompr = 1e-3;
% params.tx_ompr = linspace(1, 10, 20)'*1e-4;
ops = odeset('AbsTol', 1e-20, 'RelTol', 1e-7);
asp_grid = [.5, 1, 2, 5, 10]; %[logspace(-1, log10(2),8), linspace(2.1, 10, 7)];
params.Atot = asp_grid*1e+3;%logspace(-2, 1, 10)*1e+3;
n_t = length(params.Atot);


x0_rr = zeros(7*n_t,1);
x0_tc = zeros(6*n_t,1);

k_a_max = 0.02;%0.1;
K_d_a = 2e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr = 5;

params.tlat_ompr = params.delta*6;
params.tlat_tazc = params.tx_gfp;
params.P = 1;
T  = 1e+4;

taz_tc_ss      = zeros(length(kt_grid), length(kp_grid), n_t);
tazp_tc_ss     = zeros(length(kt_grid), length(kp_grid), n_t);
tazsum_tc_ss   = zeros(length(kt_grid), length(kp_grid), n_t);
ompr_tc_ss     = zeros(length(kt_grid), length(kp_grid), n_t);
omprp_tc_ss    = zeros(length(kt_grid), length(kp_grid), n_t);
omprsum_tc_ss  = zeros(length(kt_grid), length(kp_grid), n_t);
omprc_tc_ss    = zeros(length(kt_grid), length(kp_grid), n_t);
omprcp_tc_ss   = zeros(length(kt_grid), length(kp_grid), n_t);
omprcsum_tc_ss = zeros(length(kt_grid), length(kp_grid), n_t);
x_prot_tc_ss   = zeros(length(kt_grid), length(kp_grid), n_t);

for ii = 1 : length(kt_grid)
    for jj = 1 : length(kp_grid)
        params.kt  = kt_grid(ii);
        params.kp  = kp_grid(jj);
        params.ktc = params.kt;
        params.kpc = params.kp;
        [t_tc, x_tc] = ode15s(@(t,x)two_comp_tazcl(t,x,params),[0,T],x0_tc, ops);
        taz_tc_ss(ii,jj,:)      = x_tc(end,[0*n_t + 1 : 1*n_t]);
        tazp_tc_ss(ii,jj,:)     = x_tc(end,[1*n_t + 1 : 2*n_t]);
        tazsum_tc_ss(ii,jj,:)   = taz_tc_ss(ii,jj,:) + tazp_tc_ss(ii,jj,:);
        ompr_tc_ss(ii,jj,:)     = x_tc(end,[2*n_t + 1 : 3*n_t]);
        omprp_tc_ss(ii,jj,:)    = x_tc(end,[3*n_t + 1 : 4*n_t]);
        omprsum_tc_ss(ii,jj,:)  = ompr_tc_ss(ii,jj,:) + omprp_tc_ss(ii,jj,:);
        tazc_tc_ss(ii,jj,:)     = x_tc(end,[4*n_t + 1 : 5*n_t]);
        x_prot_tc_ss(ii,jj,:)   = x_tc(end,[5*n_t + 1 : 6*n_t]);
    end
end 

syms Taz Tazsum OmpRP OmpRsum TazC kap P tlat_ompr beta_r Asp kt kp kpc

num   = (OmpRP./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)));
tl_t  =  params.tx_gfp*num./(num + 1);
fp_tc  = [params.tlat_Taz-params.delta*Taz-kap*Taz+(Tazsum - Taz)*(kt*(OmpRsum - OmpRP)); ...
      -params.delta*OmpRP+kt*(Tazsum - Taz)*(OmpRsum - OmpRP)-(kp*Taz + kpc*TazC).*OmpRP;...
      beta_r - params.delta*OmpRsum]; 
fk_tc  = [tl_t - params.delta*TazC];
  
Ap_tc = jacobian(fp_tc, [Taz; OmpRP; OmpRsum]);
Bp_tc = jacobian(fp_tc, [TazC; beta_r]);
Bptc_a = [-k_a_max*(Asp)*Taz/(K_d_a+(Asp)); 0; 0];
Cp_tc = [0 1 0];

Ak_tc = eval(jacobian(fk_tc, [TazC]));
Bk_tc = jacobian(fk_tc, [OmpRP]);
Ck_tc = [1];

sens_fun_tc     = zeros(length(kt_grid), length(kp_grid), n_t);
r_resp_tc       = zeros(length(kt_grid), length(kp_grid), n_t);
sens_fun0_tc    = zeros(length(kt_grid), length(kp_grid), n_t);
G1_tc = cell(length(kt_grid), length(kp_grid), n_t);
G2_tc = cell(length(kt_grid), length(kp_grid), n_t);
G3_tc = cell(length(kt_grid), length(kp_grid), n_t);
K_tc  = cell(length(kt_grid), length(kp_grid), n_t);

for ii = 1 :1 :length(kt_grid)
    for jj = 1 :  length(kp_grid)
        for kk = 1 : length(params.Atot)
            Aptc_val = eval(subs(Ap_tc, [Taz; OmpRP; Tazsum; OmpRsum; TazC; kap; kt; kp; kpc], ...
                                   [taz_tc_ss(ii,jj,kk); omprp_tc_ss(ii,jj,kk); ...
                                    tazsum_tc_ss(ii,jj,kk); omprsum_tc_ss(ii,jj,kk); ...
                                    tazc_tc_ss(ii,jj,kk); params.kap_taz(kk); ...
                                    kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
            Bptc_val = eval(subs(Bp_tc, [OmpRP; kt; kp; kpc], ...
                                   [omprp_tc_ss(ii,jj,kk); ...
                                    kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
            Bptca_val = eval(subs(Bptc_a, [Taz; Asp;  kt; kp; kpc], ...
                                   [taz_tc_ss(ii,jj,kk); params.Atot(kk); ...
                                    kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
            Bktc_val = eval(subs(Bk_tc, [OmpRP; kt; kp; kpc], ...
                                   [omprp_tc_ss(ii,jj,kk); ...
                                    kt_grid(ii); kp_grid(jj); kp_grid(jj)]));
            G1_tc{ii,jj,kk} = minreal(ss(Aptc_val, Bptc_val(:,1), Cp_tc, 0));
            G2_tc{ii,jj,kk} = minreal(ss(Aptc_val, Bptc_val(:,2), Cp_tc, 0));
            K_tc{ii,jj,kk}  = ss(Ak_tc, Bktc_val, Ck_tc, 0);
            sens_fun_tc(ii,jj,kk) = norm(inv(eye(1)-K_tc{ii,jj,kk}*G1_tc{ii,jj,kk}),'inf');
            r_resp_tc(ii,jj,kk) = norm(inv(eye(1)-G1_tc{ii,jj,kk}*K_tc{ii,jj,kk})*G2_tc{ii,jj,kk},'inf');      
        end
    end
end
t1 = tic;
colors_ = {'b', 'g', 'm', 'k', 'r'}
style_ = {'',':','--','-.', ':'};
marker_ = {'', 'd', 'x', 'o'}
k_grid = 1 : length(params.Atot);
w = [0.0010,    0.0012,    0.0014,    0.0016,    0.0019,    0.0022,    0.0025,    0.0030,    0.0035,    0.0040,    0.0047,    0.0055,    0.0064,    0.0075,    0.0087,    0.0100,    0.0101,    0.0118,    0.0138,    0.0161,    0.0188,    0.0219,    0.0255,    0.0298,    0.0347,    0.0405,    0.0473,    0.0551,    0.0643,    0.0750,    0.0874,    0.1000,    0.1020,    0.1189,    0.1387,    0.1618,    0.1887,    0.2201,    0.2566,    0.2993,    0.3491,    0.4071,    0.4749,    0.5538,    0.6459,    0.7533,    0.8786,    1.0000,    1.0247,    1.1951,    1.3938,    1.6256,    1.8959,    2.2112,    2.5789,    3.0077,    3.5079,    4.0912,    4.7716,    5.5650,    6.4905,    7.5698,    8.8285,   10.0000,   10.2967,   12.0089,   14.0059,   16.3349,    19.0513,   22.2193,   25.9142,   30.2235,   35.2494,   41.1111,   47.9475,   55.9207,   65.2198,   76.0653,   88.7142, 100.0000];
figure(1)
clf;
hold on

for ii = 1 :1 
    for jj = 1 :  length(kp_grid)
        for kk2 = 3:3%1:length(k_grid)%
            kk = k_grid(kk2);
            [sv,~] = sigma(inv(eye(1)-K_tc{ii,jj,kk}*G1_tc{ii,jj,kk}),w);
            color = strcat(colors_{jj}, style_{ii});
            plot(w, sv, color, 'LineWidth',2);
        end
    end
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis([1e-3, 2*10^0, 0.25, 1.3])
toc(t1)
max(sens_fun_tc(:))
figure(2)
clf;
hold on
colors_ = {'b', 'g', 'm', 'k', 'r'}
style_ = {'',':','--','-.', ':'};
marker_ = {'', 'd', 'x', 'o'}
k_grid = 1 : length(params.Atot);
w = [0.0010,    0.0012,    0.0014,    0.0016,    0.0019,    0.0022,    0.0025,    0.0030,    0.0035,    0.0040,    0.0047,    0.0055,    0.0064,    0.0075,    0.0087,    0.0100,    0.0101,    0.0118,    0.0138,    0.0161,    0.0188,    0.0219,    0.0255,    0.0298,    0.0347,    0.0405,    0.0473,    0.0551,    0.0643,    0.0750,    0.0874,    0.1000,    0.1020,    0.1189,    0.1387,    0.1618,    0.1887,    0.2201,    0.2566,    0.2993,    0.3491,    0.4071,    0.4749,    0.5538,    0.6459,    0.7533,    0.8786,    1.0000,    1.0247,    1.1951,    1.3938,    1.6256,    1.8959,    2.2112,    2.5789,    3.0077,    3.5079,    4.0912,    4.7716,    5.5650,    6.4905,    7.5698,    8.8285,   10.0000,   10.2967,   12.0089,   14.0059,   16.3349,    19.0513,   22.2193,   25.9142,   30.2235,   35.2494,   41.1111,   47.9475,   55.9207,   65.2198,   76.0653,   88.7142, 100.0000];

for ii = 1 :1 :length(kt_grid)
    for jj = 1 :  length(kp_grid)
        for kk2 = 1:length(k_grid)
            kk = k_grid(kk2);
            [sv,~] = sigma([0 1]*inv(eye(2)-G1_tc{ii,jj,kk}*K_tc{ii,jj,kk})*G2_tc{ii,jj,kk}/d_protein,w);
            color = strcat(colors_{kk2}, style_{ii});
            plot(w, sv, color, 'LineWidth',2);
            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
        end
    end
end
toc(t1)
