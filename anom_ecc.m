function [E] = anom_ecc(M,e)

    format long

    E = M;
    k = 1;
    err = 1e-10;
    while (k>err)
        y = e*sin(E)+M;
        k = abs(abs(E)-abs(y));
        E = y;
    end
end