clear all
params.volume = 1;%602;

params.kt = 6.12*1e+3;
params.kp = 4*1e-3;
params.ktc = params.kt;
params.kpc = params.kp;

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

k_a_max = 0.02;%2;%0.1;
K_d_a = 2e+3%1e+3;
params.kap_taz = k_a_max.*(params.Atot).'./(K_d_a+(params.Atot).');

n_tr = 5;

input_ompr = 6;
input_p    = [0.1, 0.4, 0.7, 1];%linspace(0, 1, n_trc);
T  = 1e+4;
taz_rr_ss    = zeros(length(input_ompr), length(input_p), n_t);
tazp_rr_ss    = zeros(length(input_ompr), length(input_p), n_t);
tazsum_rr_ss    = zeros(length(input_ompr), length(input_p), n_t);
ompr_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprp_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprsum_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprc_rr_ss  = zeros(length(input_ompr), length(input_p), n_t);
omprcp_rr_ss = zeros(length(input_ompr), length(input_p), n_t);
omprcsum_rr_ss = zeros(length(input_ompr), length(input_p), n_t);
x_prot_rr_ss = zeros(length(input_ompr), length(input_p), n_t);
for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        params.tlat_ompr = params.delta*input_ompr(ii);
        params.P = input_p(jj);        
        [t_hill, x_hill] = ode15s(@(t,x)two_comp_sink_hill(t,x,params),[0,T],x0_rr, ops);
        taz_rr_ss(ii,jj,:)      = x_hill(end,[0*n_t + 1 : 1*n_t]);
        tazp_rr_ss(ii,jj,:)     = x_hill(end,[1*n_t + 1 : 2*n_t]);
        tazsum_rr_ss(ii,jj,:)   = taz_rr_ss(ii,jj,:) + tazp_rr_ss(ii,jj,:);
        ompr_rr_ss(ii,jj,:)     = x_hill(end,[2*n_t + 1 : 3*n_t]);
        omprp_rr_ss(ii,jj,:)    = x_hill(end,[3*n_t + 1 : 4*n_t]);
        omprsum_rr_ss(ii,jj,:)  = ompr_rr_ss(ii,jj,:) + omprp_rr_ss(ii,jj,:);
        omprc_rr_ss(ii,jj,:)    = x_hill(end,[4*n_t + 1 : 5*n_t]);
        omprcp_rr_ss(ii,jj,:)   = x_hill(end,[5*n_t + 1 : 6*n_t]);
        omprcsum_rr_ss(ii,jj,:) = omprc_rr_ss(ii,jj,:) + omprcp_rr_ss(ii,jj,:);
        x_prot_rr_ss(ii,jj,:)   = x_hill(end,[6*n_t + 1 : 7*n_t]);
    end
end


syms X Taz Tazsum deltaTazP deltaTaz deltaOut TazC OmpRP OmpRsum OmpRC OmpRCsum kap P tlat_ompr beta_r Asp w

num   = ((OmpRP)./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)));
tl_rc  = P*tl_x;

fp_rr  = [params.tlat_Taz-params.delta*Taz-kap*Taz+(Tazsum - Taz+deltaTazP)*(params.kt*(OmpRsum - OmpRP)+params.ktc*OmpRC); ...
	      params.tlat_Taz-params.delta*Tazsum;...
         -params.delta*OmpRP+params.kt*(Tazsum - Taz+deltaTazP)*(OmpRsum - OmpRP)-params.kp*(Taz+deltaTaz).*OmpRP;...
          beta_r - params.delta*OmpRsum]; 
fk_rr  = [tl_rc+deltaOut-params.delta*OmpRC-params.ktc*(Tazsum - Taz+deltaTazP)*OmpRC+params.kpc*(Taz+deltaTaz)*(OmpRCsum - OmpRC);... 
          tl_rc+deltaOut-params.delta*OmpRCsum];

Acl = jacobian([fp_rr; fk_rr; tl_x+deltaOut-params.delta*X], [Taz; Tazsum; OmpRP; OmpRsum; OmpRC; OmpRCsum; X]);
Bcl = jacobian([fp_rr; fk_rr; tl_x+deltaOut-params.delta*X], [deltaOut; deltaTaz; deltaTazP]);
Ccl = [0 0 0 0 0 0 1];


sens_fun_rr     = zeros(length(input_ompr), length(input_p), n_t);
Scl = cell(length(input_ompr), length(input_p), n_t);



for ii = 1 : length(input_ompr)
    for jj = 1 : length(input_p)
        for kk = 1 : length(params.Atot)
            Acl_val = eval(subs(Acl, [Taz; OmpRP; Tazsum; OmpRsum; OmpRC;  OmpRCsum; X; kap; P; deltaOut; deltaTazP; deltaTaz], ...
                                   [taz_rr_ss(ii,jj,kk); omprp_rr_ss(ii,jj,kk); ...
                                    tazsum_rr_ss(ii,jj,kk); omprsum_rr_ss(ii,jj,kk); ...
                                    omprc_rr_ss(ii,jj,kk); omprcsum_rr_ss(ii,jj,kk); x_prot_rr_ss(ii,jj,kk); ...
                                    params.kap_taz(kk); input_p(jj); 0; 0; 0]));
            Bcl_val = eval(subs(Bcl, [Taz; OmpRP; Tazsum; OmpRsum; OmpRC; OmpRCsum; X; kap; P; deltaOut; deltaTazP; deltaTaz], ...
                                   [taz_rr_ss(ii,jj,kk); omprp_rr_ss(ii,jj,kk); ...
                                    tazsum_rr_ss(ii,jj,kk); omprsum_rr_ss(ii,jj,kk); ...
                                    omprc_rr_ss(ii,jj,kk); omprcsum_rr_ss(ii,jj,kk); x_prot_rr_ss(ii,jj,kk); ...
                                    params.kap_taz(kk); input_p(jj); 0; 0; 0]));
            Scl{ii,jj,kk} = minreal(ss(Acl_val, Bcl_val, Ccl, 0));
            sens_fun_rr(ii,jj,kk)    = norm(Scl{ii,jj,kk},'inf');
%             sens_fun0_rr(ii,jj,kk)   = freqresp(inv(eye(1)-K_rr{ii,jj,kk}*G1_rr{ii,jj,kk}),0);
%             r_resp_rr(ii,jj,kk)      = norm([0 1]*inv(eye(2)-G1_rr{ii,jj,kk}*K_rr{ii,jj,kk})*G2_rr{ii,jj,kk},'inf');
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

for ii = 1 :1 % length(input_ompr)
    for jj = 1 :  length(input_p)
        for kk2 = 3:3%length(k_grid):length(k_grid)
            kk = k_grid(kk2);
%     	    sigma(Scl{ii,jj,kk});
%             keyboard
            [sv,~] = sigma(Scl{ii,jj,kk}(1),w);
            color = strcat(colors_{jj}, style_{ii});
            plot(w, sv(1,:), color, 'LineWidth',2);
            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
        end 
    end
end
% axis([1e-3, 10^2, 0.25, 1.3])
max(sens_fun_rr(:))
toc(t1)
