function localExt = SIRS_solveOnce_ie_op(x0, t0, tf, k, event_type)

x            = x0;
t            = t0;
localExt     = x(2);
stoch_matrix = [-1 1 0; 0 -1 1; 1 0 -1];

while(t < tf)
    % compute propensity
    b    = k.*x;
    b(1) = b(1) * x(2);
    b0   = sum(b);
    
	% time to the next reaction
	tau = log(1/rand) / b0;
	t   = t + tau;
	
	if(t <= tf)
        j = find(cumsum(b)/b0 > rand, 1);
		x = x + stoch_matrix(j,:);
				
		currentEevntVal = x(2);
		% record extreme values
		if (event_type == 1)
			localExt = max(localExt, currentEevntVal);
		else
			localExt = min(localExt, currentEevntVal);
		end
	end
end