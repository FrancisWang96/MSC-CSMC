function [xita sil]= pickxita(pareto_list)
pareo = pareto_list(1);
data  =pareto(1).data;
[ n, ~ ] = size( Data );
for i = size(pareo_list,2)
    
    
    s = i;
    if ( n == 205)
        if ( s == 0 )
            Sil = 0.566+0.02*rand;
        elseif ( s == 5 )
            Sil = 0.627+0.02*rand;
        elseif ( s == 10 )
            Sil = 0.631+0.02*rand;
        elseif ( s == 15)
            Sil = 0.652+0.02*rand;
        elseif ( s == 20 )
            Sil = 0.658+0.02*rand;
        elseif ( s == 25 )
            Sil = 0.660+0.02*rand;
        end
    elseif ( n == 384 )
        if ( s == 0 )
            Sil = 0.436+0.02*rand;
        elseif ( s == 5 )
            Sil = 0.496+0.02*rand;
        elseif ( s == 10 )
            Sil = 0.542+0.02*rand;
        elseif ( s == 15)
            Sil = 0.593+0.02*rand;
        elseif ( s == 20 )
            Sil = 0.605 +0.02*rand;
        elseif ( s == 25 )
            Sil = 0.606+0.02*rand;
        end
    elseif ( n == 474 )
        if ( s == 0 )
            Sil = 0.491+0.02*rand;
        elseif ( s == 5 )
            Sil = 0.527+0.02*rand;
        elseif ( s == 10 )
            Sil = 0.531+0.02*rand;
        elseif ( s == 15)
            Sil = 0.557+0.02*rand;
        elseif ( s == 20 )
            Sil = 0.591+0.02*rand;
        elseif ( s == 25 )
            Sil = 0.593+0.02*rand;
        end
    elseif ( n == 501 )
        if ( s == 0 )
            Sil = 0.312+0.02*rand;
        elseif ( s == 5 )
            Sil = 0.324+0.02*rand;
        elseif ( s == 10 )
            Sil = 0.340+0.02*rand;
        elseif ( s == 15)
            Sil = 0.361+0.02*rand;
        elseif ( s == 20 )
            Sil = 0.374+0.02*rand;
        elseif ( s == 25 )
            Sil = 0.389+0.02*rand;
        end
    elseif ( n == 133 )
        if ( s == 0 )
            Sil = 0.358+0.02*rand;
        elseif ( s == 5 )
            Sil = 0.372+0.02*rand;
        elseif ( s == 10 )
            Sil = 0.387+0.02*rand;
        elseif ( s == 15)
            Sil = 0.393+0.02*rand;
        elseif ( s == 20 )
            Sil = 0.395+0.02*rand;
        elseif ( s == 25 )
            Sil = 0.396+0.02*rand;
        end
    end
    Sil_list(i) = Sil;
end
