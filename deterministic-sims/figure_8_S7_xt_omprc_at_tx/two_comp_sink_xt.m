function [dx] = two_comp_sink_xt(t, x, params, dist_signals)

if (size(x,1) == 1)
    x = x(:);
end
n_sim = size(x,1)/9;

Taz         = x(0*n_sim + 1  :  1*n_sim);
TazP        = x(1*n_sim + 1  :  2*n_sim);
OmpR        = x(2*n_sim + 1  :  3*n_sim);
OmpRP       = x(3*n_sim + 1  :  4*n_sim);
OmpRC       = x(4*n_sim + 1  :  5*n_sim);
OmpRCP      = x(5*n_sim + 1  :  6*n_sim);
X           = x(6*n_sim + 1  :  7*n_sim);
Ph          = x(7*n_sim + 1  :  8*n_sim);
HK          = x(8*n_sim + 1  :  9*n_sim);

num   = (OmpRP./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
dist_signals_t = dist_signals(t);
tl_x  = (params.tx_gfp*(num./(num + 1)))+dist_signals_t(4);
tl_r = params.tlat_ompr+dist_signals_t(5);
tl_t = params.tlat_Taz + dist_signals_t(6);
tl_rc  = params.P*tl_x; %params.tlat_omprc*num./(num + 1);
kup = params.kup_hill+dist_signals_t(1);
dx = [tl_t-params.delta_Taz*Taz-kup.*Taz+TazP.*(params.kt*OmpR+params.ktc*OmpRC+dist_signals_t(2)); ...
      kup.*Taz-TazP.*(params.kt*OmpR+params.ktc*OmpRC+dist_signals_t(2))-params.delta_Taz*TazP;...
      %
      tl_r-params.delta_OmpR*OmpR-(params.kt.*TazP+params.ktx.*HK+dist_signals_t(2)).*OmpR+(params.kp.*Taz+params.kpx.*Ph+dist_signals_t(3)).*OmpRP;... -
      -params.delta_OmpR*OmpRP+(params.kt.*TazP+params.ktx.*HK+dist_signals_t(2)).*OmpR-(params.kp.*Taz+params.kpx.*Ph+dist_signals_t(3)).*OmpRP-params.k_titr*OmpRCP.*OmpRP;...%;... 
      %
      tl_rc-params.delta_OmpR*OmpRC-(params.ktc.*TazP+params.ktx.*HK+dist_signals_t(2)).*OmpRC+(params.kpc.*Taz+params.kpx.*Ph+dist_signals_t(3)).*OmpRCP;... -
      -params.delta_OmpR*OmpRCP+(params.ktc.*TazP+params.ktx.*HK+dist_signals_t(2)).*OmpRC-(params.kpc.*Taz+params.kpx.*Ph+dist_signals_t(3)).*OmpRCP-params.k_titr*OmpRCP.*OmpRP;...%; ...
      %
      tl_x - params.delta_OmpR*X;...
      params.tlat_Taz-params.delta_Taz*Ph-params.kup_xt.*Ph+HK.*(params.ktx*OmpR+params.ktx*OmpRC); ...
      params.kup_xt.*Ph-HK.*(params.ktx*OmpR+params.ktx*OmpRC)-params.delta_Taz*HK];
end
