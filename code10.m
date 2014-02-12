% Vignesh Code Adpated Version
M = size(E,1);
N = size(E,2);

En = E_noise;

p_old = p_init;
q_old = q_init;
p_new = zeros(size(En));
q_new = zeros(size(En));

maxiter = 300; %Maximmum number of iterations 


isize= M;
iter = 0;
lambda =   50;

while (1)       
    %For p
    if(1)
        for i = 2:isize-1
            for j = 2:isize-1
                if mask(i,j) == 0
                    continue;
                end
                if boundary(i,j) == 0     

                    p_new(i,j) = Ravg(p_old, i, j)  + (1/lambda)*( En(i,j) - Rval(p_old(i,j),q_old(i,j),s) )*Rp(p_old(i,j),q_old(i,j),s(1),s(2)) ;                   
                else
                    p_new(i,j) = p_old(i,j);
                end            
            end
        end
    end
    %For q
    if(1)

         for i = 2:isize-1
            for j = 2:isize-1
                if mask(i,j) == 0
                    continue;
                end
                if boundary(i,j) == 0
                    q_new(i,j) = Ravg(q_old, i, j)  + (1/lambda)*( En(i,j) - Rval(p_old(i,j),q_old(i,j),s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2)) ;            
                else
                    q_new(i,j) = q_old(i,j);
                end
            end
         end
    end
    
    if  (iter == maxiter)
        break;
        
    else
        disp(iter)
    end

    p_old = p_new;
    q_old = q_new;
    
    iter = iter + 1;
end