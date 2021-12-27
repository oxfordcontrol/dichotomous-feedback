function [dx] = two_comp_sink_xt(t, x, params, w)

if (size(x,1) == 1)
    x = x(:);
end
n_sim = size(x,1)/7;

Taz         = x(0*n_sim + 1  :  1*n_sim);
TazP        = x(1*n_sim + 1  :  2*n_sim);
OmpR        = x(2*n_sim + 1  :  3*n_sim);
OmpRP       = x(3*n_sim + 1  :  4*n_sim);
OmpRC       = x(4*n_sim + 1  :  5*n_sim);
OmpRCP      = x(5*n_sim + 1  :  6*n_sim);
X           = x(6*n_sim + 1  :  7*n_sim);

% dist_signals_t = ;
dist_signal = 1e-6*sin(w*t);
num   = ((OmpRP+dist_signal)./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)));
tl_r = params.tlat_ompr;
tl_t = params.tlat_Taz;
tl_rc  = params.P*tl_x; %params.tlat_omprc*num./(num + 1);
kup = params.kap_taz;
dx = [tl_t-params.delta*Taz-kup.*Taz+TazP.*(params.kt*OmpR+params.ktc*OmpRC); ...
      kup.*Taz-TazP.*(params.kt*OmpR+params.ktc*OmpRC)-params.delta*TazP;...
      %
      tl_r-params.delta*OmpR-(params.kt.*TazP).*OmpR+(params.kp.*Taz).*OmpRP;... -
      -params.delta*OmpRP+(params.kt.*TazP).*OmpR-(params.kp.*Taz).*OmpRP;...%;... 
      %
      tl_rc-params.delta*OmpRC-(params.ktc.*(TazP)).*OmpRC+(params.kpc.*(Taz)).*OmpRCP;... -
      -params.delta*OmpRCP+(params.ktc.*(TazP)).*OmpRC-(params.kpc.*(Taz)).*OmpRCP;...%; ...
      %
      tl_x - params.delta*X];
end
