function localExt = RI_solveOnce_ie_op(x0, t0, tf, k, event_type)

x            = x0;
t            = t0;
localExt     = x(2);

while(t < tf)    
    % compute propensity
    b  = k.*x;
    b0 = b(1) + b(2);
      
    % time to the next reaction
    tau = log(1/rand) / b0;
	t   = t + tau;
	
	if(t <= tf)
        if cumsum(b)/b0 > rand
            x = x + [-1 1];
        else
            x = x + [1 -1];
        end
		
		% check for the target event
		currentEevntVal = x(2);
		% record extreme values
		if (event_type == 1)
			localExt = max(localExt, currentEevntVal);
		else
			localExt = min(localExt, currentEevntVal);
		end
	end
end