function [ u,t ] = generateSquareWave( tau, Tf, Ts, numOfWaves)
      
        [ u,t ] = gensig('square',tau, Tf, Ts);
        for i=1:length(u)/2
          u(i) = 1;
        end

        t_new = 0:10^-5:0.7*10^-3;
        u_new = zeros(length(t_new),1);
        for i=1:length(t_new)

           if(i<length(t))
               u_new(i) = u(i);
           else
               u_new(i) = 0;
           end
           
        end
        
        u = u_new;
        t = t_new;
        
        t_new = zeros(numOfWaves*length(t),1);
        u_new = zeros(length(t_new),1);
        
        step = 0;
        offset = 0;
        for i=1:numOfWaves
            for j=1:length(t)
               u_new(j+offset) = u(j);
               t_new(j+offset) = step;
               step = step+10^-5;
            end
            offset = offset+length(t);
        end
        
        u = u_new;
        t = t_new;

end

