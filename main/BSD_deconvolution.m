function [N,C,CostFunction]  = BSD_deconvolution(Fluorescence,Ninit,b...
                ,gamma,delta,eta,lambda,mu,conditioner, slowOptimizer, superResolution, I,I1,M,d0,d1,d2)
        % This function solves the problem:
        % N=argmin_{N' \in R_+} sum( 0.5*(F-K*N'-b)^2+lambda*N' ).
        % Inputs: Fluorescence: a single neuron fluorescence signal.
        %         params: The usual parameters. We use here
        %         gamma,delta,eta,b (for the squared error), lambda (for
        %         the penalty).
        %         p: The neuron number.
        %         I,M,d0,d1,d2: I & M are initialized sparse matrices. d0,
        %         d1, d2 are the target indices.
        % Output: n: the inferred spikes (in inferred units)
        %         C: the corresponding fluorescence trace (obtained by
        %         convolution of n, not normalized, no offset).
        H1 = I; % initialize memory for hessian.
        M(d0) = eta;
        M(d1) =  -(gamma+delta) * eta; 
        M(d2) =  delta * eta; % The matrix such that n=M*C (The discrete Differential Equation).
        
        N=0*Ninit+0.01;
        C= M\N;
       z = max(1,sqrt(max(lambda)));             % weight on barrier function
       if slowOptimizer % Slow but more accurate optimizer, usefull sometimes for SR
           zmin = 1e-15;
           zreason = 2;
       else
           zmin = 1e-10;
           zreason = 10;
       end;
       while z>zmin                        % this is an arbitrary threshold
            Error=Fluorescence-C(superResolution:superResolution:end)-b;
            CostFunction = 0.5*sum(Error.^2)+ 0.5 * conditioner * sum(mu .* N.^2)+ sum(lambda .* N)-z*sum(log(N));
            s    = 1;                           % step size
            d    = Inf;                           % direction
            count = 0;
             while (norm(d)>1e-2 && s > 1e-2) && ( superResolution==1 || count<100)      % converge for this z (again, these thresholds are arbitrary)
                 count = count+1;
                % Compute Gradient
                Gradient    = M'*(lambda+ conditioner * (mu .*N) -z * N.^-1);
                Gradient(superResolution:superResolution:end) = Gradient(superResolution:superResolution:end)-Error;             
                
                
                % Compute Hessian
                H1(d0)  = z*N.^-2+conditioner*mu;                % log barrier part of the Hessian
                H       = I1 + (M'*H1*M);     % Hessian                
                % Compute Direction and step.
                d   = -H\Gradient;                     % direction to step using newton-raphson
                hit = -N./(M*d);                % step within constraint boundaries
                hit(hit<0)=[];                  % ignore negative hits
                if any(hit<1)
                    s = min(1,0.99*min(hit(hit>0)));
                else
                    s = 1;
                end
                CostFunctionLoop = CostFunction + abs(CostFunction);
                CLoop =C;
                while CostFunctionLoop>=CostFunction + abs(CostFunction)*(1e-4)         % make sure newton step doesn't increase objective
                    CLoop  = C+s*d;
                    N   = M*CLoop;
                    Error=Fluorescence-CLoop(superResolution:superResolution:end)-b;
                    CostFunctionLoop = 0.5*sum(Error.^2) + 0.5 * conditioner * sum(mu .* N.^2)+ sum(lambda .* N) -z*sum(log(N));
                    s   = s/5;                  % if step increases objective function, decrease step size
                    if s<1e-5; break; end      % if decreasing step size just doesn't do it
                end
                C    = CLoop;                      % update C
                CostFunction = CostFunctionLoop;                   % update post
             end
            if (superResolution >1 && count>=100)
                if conditioner>0
                    conditioner = conditioner * 10;
                else
                    conditioner = 1e-6;
                end;
                display(sprintf('Newton did not converge for z=%.e ;increase to conditioner=%.e and reset z.',z,conditioner));
                z = max(1,sqrt(max(lambda)));
                N=0*Ninit+0.01;  
                C= M\N;                                
            elseif max(isnan(N));
                zmin = z*zreason;
                display(sprintf('Found NaN for z = %.e; increase zmin to %.e and reset z.',z,zmin));
                z = max(1,sqrt(max(lambda)));
                N=0*Ninit+0.01;
                C= M\N;                
            else
                z=z/zreason;                             % reduce z (sequence of z reductions is arbitrary)
            end;
       end
       Error = Fluorescence - C(superResolution:superResolution:end)-b;
       CostFunction = 0.5*sum(Error.^2) + 0.5 * conditioner * sum(mu .* N.^2)+ sum(lambda .* N);
    end






 