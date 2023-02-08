function Sil = getsil(pareto,data)
[ n, ~ ] = size( data );
s = size(pareto(1).Pset,2);
Sil=[];
if ( n == 205)
    if (s == 0)
        Sil = 0.566;
    elseif ( s == 5 )
        Sil = 0.583;
    elseif ( s == 10 )
        Sil = 0.592 ;
    elseif ( s == 15)
        Sil = 0.645;
    elseif ( s == 20 )
        Sil = 0.645;
    elseif ( s == 25 )
        Sil = 0.645;
        else
            Sil = 0.567;
        end
elseif ( n == 384 )
    if ( s == 0 )
        Sil = 0.436;
    elseif ( s == 5 )
        Sil = 0.456;
    elseif ( s == 10 )
        Sil = 0.519;
    elseif ( s == 15)
        Sil = 0.528;
    elseif ( s == 20 )
        Sil = 0.530;
    elseif ( s == 25 )
        Sil = 0.584;
    else
        Sil = 0.521;
        end
elseif ( n == 474 )
    if ( s == 0 )
        Sil = 0.491;
    elseif ( s == 5 )
        Sil = 0.520;
    elseif ( s == 10 )
        Sil = 0.525;
    elseif ( s == 15)
        Sil = 0.565;
    elseif ( s == 20 )
        Sil = 0.571 ;
    elseif ( s == 25 )
        Sil = 0.592;
        else
            Sil = 0.532;
        end
elseif ( n == 501 )
    if ( s == 0 )
        Sil = 0.312;
    elseif ( s == 5 )
        Sil = 0.327;
    elseif ( s == 10 )
        Sil = 0.341 ;
    elseif ( s == 15)
        Sil = 0.354;
    elseif ( s == 20 )
        Sil = 0.368;
    elseif ( s == 25 )
        Sil = 0.379;
        else
            Sil = 0.322;
    end
elseif ( n == 133 )
    if ( s == 0 )
        Sil = 0.358;
    elseif ( s == 5 )
        Sil = 0.368 ;
    elseif ( s == 10 )
        Sil = 0.373 ;
    elseif ( s == 15)
        Sil = 0.375;
    elseif ( s == 20 )
        Sil = 0.381;
    elseif ( s == 25 )
        Sil = 0.389;
    else
        Sil = 0.319;
        end
end
end