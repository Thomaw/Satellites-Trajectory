function [Ascension_droite_NA_deg,Ascension_droite_NA,Arg_p,Nu,i]=...
    radian_converter(Ascension_droite_NA,Arg_p,Nu,i)
    

%-----------------------------------------------------%
%-------------- Conversions nécessaires --------------%
% Conversion des valeurs en radian
Ascension_droite_NA_deg = Ascension_droite_NA; % utile dans la suite pour greenwich0
Ascension_droite_NA = Ascension_droite_NA*pi/180;
Arg_p               = Arg_p*pi/180;
Nu                  = Nu*pi/180;
i                   = i*pi/180;
%-----------------------------------------------------%
%-----------------------------------------------------%

end