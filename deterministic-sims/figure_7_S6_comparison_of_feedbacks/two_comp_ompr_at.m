function [dx] = two_comp_tazc(t, x, params)

if (size(x,1) == 1)
    x = x(:);
end
n_sim = size(x,1)/6;

Taz         = x(0*n_sim + 1  :  1*n_sim);
TazP        = x(1*n_sim + 1  :  2*n_sim);
OmpR        = x(2*n_sim + 1  :  3*n_sim);
OmpRP       = x(3*n_sim + 1  :  4*n_sim);
FC          = x(4*n_sim + 1  :  5*n_sim);
X           = x(5*n_sim + 1  :  6*n_sim);

num   = (OmpRP./params.Kdr).^params.hill_coeff; %  params.tx_gfp*OmpRP;%
tl_x  = (params.tx_gfp*(num./(num + 1)));
tl_f = params.P*tl_x;
tl_r  = params.tlat_ompr;
dx = [params.tlat_Taz-params.delta*Taz-params.kap_taz.*Taz+params.kt*TazP.*OmpR; ...
     -params.delta*TazP+params.kap_taz.*Taz-params.kt*TazP.*OmpR;...
      tl_r-params.delta*OmpR-params.kt.*TazP.*OmpR+(params.kp.*Taz).*OmpRP-params.kpc.*FC.*OmpR;... 
      -params.delta*OmpRP+params.kt.*TazP.*OmpR-(params.kp.*Taz+params.kpc.*FC).*OmpRP;... 
      tl_f - params.delta*FC-params.kpc.*FC.*OmpRP-params.kpc.*FC.*OmpR;...
      tl_x - params.delta*X];
end
