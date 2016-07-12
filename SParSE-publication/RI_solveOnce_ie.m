function [counter_out, localExt] = RI_solveOnce_ie(x0, t0, tf, k, gammas, target_e, event_type)

counter_out = 0;
localExt    = x0(2);

x           = x0;
t           = t0;
flag        = (x(2)==target_e);

while(t < tf)
    % check for a rare event
    if flag
        counter_out = 1;
        break;
    end
    
    % compute propensity
    b  = k.*x.*gammas;
    b0 = sum(b);
      
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
		flag = (x(2)==target_e);
        
		% record extreme values
		if (event_type == 1)
			localExt = max(localExt, x(2));
		else
			localExt = min(localExt, x(2));
		end
	end
end