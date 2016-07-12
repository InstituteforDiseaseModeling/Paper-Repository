function gammai = interpolateGammas(gammas_opt, dists)

signs = sign(dists);
sprod = signs(1:(end-1)).*signs(2:end);
sind = find(sprod==-1,1,'first');

% weighted average according to the distance from the target
d_tot = sum(abs(dists(sind:sind+1)));
gammai = gammas_opt(:,sind)'*abs(dists(sind))/d_tot + gammas_opt(:,sind+1)'*abs(dists(sind+1))/d_tot;