function [dx] = two_comp_sink_hill(t, x, params)

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

num   = (OmpRP./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)));
tl_rc  = params.P*tl_x;
dx = [params.tlat_Taz-params.delta*Taz-params.kap_taz.*Taz+TazP.*(params.kt*OmpR+params.ktc*OmpRC); ...
      params.kap_taz.*Taz-TazP.*(params.kt*OmpR+params.ktc*OmpRC)-params.delta*TazP;...
      params.tlat_ompr-params.delta*OmpR-params.kt.*TazP.*OmpR+params.kp.*Taz.*OmpRP;... -
      -params.delta*OmpRP+params.kt.*TazP.*OmpR-params.kp.*Taz.*OmpRP;... 
      tl_rc-params.delta*OmpRC-params.ktc.*TazP.*OmpRC+params.kpc.*Taz.*OmpRCP;... -
      -params.delta*OmpRCP+params.ktc.*TazP.*OmpRC-params.kpc.*Taz.*OmpRCP; ...
      tl_x - params.delta*X];
end
