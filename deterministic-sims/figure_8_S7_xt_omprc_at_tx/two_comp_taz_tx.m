function [dx] = two_comp_tazc(t, x, params)

if (size(x,1) == 1)
    x = x(:);
end
n_sim = size(x,1)/8;

Taz         = x(0*n_sim + 1  :  1*n_sim);
TazP        = x(1*n_sim + 1  :  2*n_sim);
OmpR        = x(2*n_sim + 1  :  3*n_sim);
OmpRP       = x(3*n_sim + 1  :  4*n_sim);
FC          = x(4*n_sim + 1  :  5*n_sim);
Y           = x(5*n_sim + 1  :  6*n_sim);
Ph          = x(6*n_sim + 1  :  7*n_sim);
HK          = x(7*n_sim + 1  :  8*n_sim);

num   = (OmpRP./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_y  = (params.tx_gfp*(num./(num + 1)));
tl_f = params.P*tl_y;
tl_r  = params.tlat_ompr;
tl_t  = params.tlat_Taz./(1 + (FC./params.Kdfc).^1);
%keyboard
dx = [tl_t-params.delta*Taz-params.kap_taz.*Taz+params.kt*TazP.*OmpR; ...
     -params.delta*TazP+params.kap_taz.*Taz-params.kt*TazP.*OmpR;...
      tl_r-params.delta*OmpR-(params.kt.*TazP+params.ktx.*HK).*OmpR+(params.kp.*Taz+params.kpx.*Ph).*OmpRP;... 
      -params.delta*OmpRP+(params.kt.*TazP+params.ktx.*HK).*OmpR-(params.kp.*Taz+params.kpx.*Ph).*OmpRP;... 
      tl_f - params.delta*FC;...
      tl_y - params.delta*Y;...
      params.tlat_Taz-params.delta_Taz*Ph-params.kup_xt.*Ph+HK.*(params.ktx*OmpR); ...
      params.kup_xt.*Ph-HK.*(params.ktx*OmpR)-params.delta_Taz*HK];
end
