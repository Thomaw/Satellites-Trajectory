clear;
close all;
clc;


%% Caractéristique objet célèste

fprintf('Mercury \nVenus \nEarth \nMoon \nMars \nJupiter \nSaturn \nUranus \nNeptune\n\n')
ch=input('Rentrez le nom de la planète que vous souhaité : ','s');

%{ 
On définit :
R : Rayon de la Terre                     [km]
mu : Paramètre gravitationnel standard    [km^3/s^2]
Omega : Vitesse de rotation terrestre     [rad/s]
%}

switch ch
    case 'Mercury'
        image_file = 'Planets\Mercury_map.jpg';
        R_e = 2440;
        mu_e = 22030;
        Omega_e = (2*pi/(3600*24*58.65));
        
    case 'Venus'
        image_file = 'Planets\Venus_map.jpg';
        R_e = 6052;
        mu_e = 324900;
        Omega_e = -(2*pi/(3600*24*243)); % Retrograde(on met donc un moins)
        
    case 'Earth'
        image_file = 'Planets\Earth_map.jpg';
        R_e = 6378;
        mu_e = 398600;
        Omega_e = (2*pi/(3600*23.9345));
        
    case 'Moon'
        image_file = 'Planets\Moon_map.jpg';
        R_e = 1737;
        mu_e = 4903;
        Omega_e = (2*pi/(3600*24*27.32));
        
    case 'Mars'
        image_file = 'Planets\Mars_map.jpg';
        R_e = 3396;
        mu_e = 42828;
        Omega_e = (2*pi/(3600*24.62));
        
    case 'Jupiter'
        image_file = 'Planets\Jupiter_map.jpg';
        R_e = 71490;
        mu_e = 126686000;
        Omega_e = (2*pi/(3600*9.925));
        
    case 'Saturn'
        image_file = 'Planets\Saturn_map.jpg';
        R_e = 60270;
        mu_e = 37931000;
        Omega_e = (2*pi/(3600*10.66));
        
    case 'Uranus'
        image_file = 'Planets\Uranus_map.jpg';
        R_e = 25560;
        mu_e = 5794000;
        Omega_e = -(2*pi/(3600*17.24)); % Retrograde(on met donc un moins)
        
    case 'Neptune'
        image_file = 'Planets\Neptune_map.jpg';
        R_e = 24760;
        mu_e = 6835100;
        Omega_e = (2*pi/(3600*16.11));
        
    otherwise
        fprintf('Vous n''avez pas choisi de planètes. \nVeuillez réessayer\n\n')
        return
end
     



fprintf('\n----------- Caratéristiques de base ------------ \n')
fprintf('-------------------------------------------------- \n')
fprintf('Rayon de la planète = %f Km \n', R_e)
fprintf('Paramètre gravitationnel standard = %f km^3/s^2 \n', mu_e)
fprintf('Vitesse de rotation = %f rad/s \n', Omega_e)
fprintf('-------------------------------------------------- \n')
fprintf('-------------------------------------------------- \n\n')
% La rotation de la Terre (Omega_e) se situe sur l'axe Z.



%% Définition des couleurs de plot

colors={[1 0 0] [0 1 0] [1 1 0] [1 0 1] [1 0.27 0] [0 1 1] [0.33 0 0.51]...
    [0.545 0.27 0.075] [0.5 0 0] [0.86 0.08 0.235] [0 0 0.5] [0.5 0.5 0]...
    [1 0.98 0.8] [0.82 0.4 0.12] [1 0.41 0.71] [0.12 0.56 1] [0.4 0.8 0.67]...
    [0.6 0.98 0.6] [0.75 0.75 0.75] [1 0.08 0.56]};



%% Demande du nombre de satellites

k = -1;
fprintf('Rentrez le nombre de satellites :\n')
while k<=0 || floor(k)~=k
    k=input('');
end


%% Initialisation des paramètres

Ascension_droite_NA_deg=cell(1,k);
Ascension_droite_NA=cell(1,k);
Arg_p=cell(1,k);
Nu=cell(1,k);
i=cell(1,k);
a=cell(1,k);
e=cell(1,k);

rp=cell(1,k);
ra=cell(1,k);
Vp=cell(1,k);
Va=cell(1,k);
n=cell(1,k);
p=cell(1,k);
T=cell(1,k);
h=cell(1,k);
h1=cell(1,k);
h2=cell(1,k);
h3=cell(1,k);
n1=cell(1,k);
n2=cell(1,k);
n3=cell(1,k);
N=cell(1,k);

heures=cell(1,k);
minutes=cell(1,k);
secondes=cell(1,k);

norb=cell(1,k);
delta_T=cell(1,k);

tp=cell(1,k);
Nu_inf=cell(1,k);
V_inf=cell(1,k);

sinF0=cell(1,k);
cosF0=cell(1,k);
F0=cell(1,k);

cos_E0=cell(1,k);
sin_E0=cell(1,k);
E0=cell(1,k);

mmm=cell(1,k);

ellipse_per=cell(1,k);

%% Boucle calcul

for aa=1:k
    fprintf('\nSelectionnez l''orbite de l''objet %d que vous souhaitez \n',aa)
    fprintf('1) Elliptique \n2) Parabolique \n3) Hyperbolique\n')
    
    mmm{aa}=0;
    while (mmm{aa}~=1 && mmm{aa}~=2 && mmm{aa}~=3)
        mmm{aa}=input('');
    end
    
    
    if mmm{aa}==1
        %trajectoire elliptique
        
        fprintf('\nRentrez l''Ascension droite du nœud ascendant entre [0,360[ en degré:\n')
        Ascension_droite_NA{aa} = -1;
        while Ascension_droite_NA{aa}<0 || Ascension_droite_NA{aa}>360 
            Ascension_droite_NA{aa} = input('');
        end
        
        fprintf('Rentrez l''argument du périgée entre [0,360[ en degré:\n')
        Arg_p{aa} = -1;
        while Arg_p{aa}<0 || Arg_p{aa}>360
            Arg_p{aa}=input('');
        end
        
        fprintf('Rentrez l''anomalie vraie (pour t=0) entre [0,360[ en degré:\n')
        Nu{aa} = -1;
        while Nu{aa}<0 || Nu{aa}>360
            Nu{aa} = input('');
        end
        
        fprintf('Rentrez l''inclinaison entre [-90, 90] en degré:\n')
        i{aa} = -91;
        while i{aa}<-90 || i{aa}>90
            i{aa} = input('');
        end
        
        fprintf('Rentrez le demi-grand axe de la trajectoire (>%f) en Km:\n',R_e)
        a{aa} = 0;
        while a{aa}<R_e
            a{aa} = input('');
        end
        
        ecc_max = 1-R_e/a{aa};
        fprintf('Rentrez l''excentricité. Il ne faut pas dépasser %f :\n', ecc_max)
        e{aa} = -1;
        while e{aa}>ecc_max || e{aa}<0
            e{aa} = input('');
        end
        
        
        
        
        fprintf('\n\n---------- Caracteristiques de l''orbite %d--------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        
        fprintf('Ascension droite du nœud ascendant numéro %d = %f deg \n',aa,Ascension_droite_NA{aa})
        fprintf('Argument au périgée numéro %d = %f deg \n',aa,Arg_p{aa})
        fprintf('Anomalie vraie au départ numéro %d = %f deg \n',aa,Nu{aa})
        fprintf('Inclinaison numéro %d = %f deg \n',aa,i{aa})
        fprintf('Demi-grand axe numéro %d = %f en Km \n',aa,a{aa})
        fprintf('Excentricité numéro %d = %f deg \n',aa,e{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n')
        
        
        
        
        [Ascension_droite_NA_deg{aa},Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa}]=...
            radian_converter(Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa});
        
        
        
        %% Calcul des paramètres orbitaux
        
        rp{aa} = a{aa}*(1-e{aa});                % rayon périgée [km]
        ra{aa} = a{aa}*(1+e{aa});                % rayon apogée  [km]
        
        Vp{aa} = sqrt(mu_e*(2/rp{aa}-1/a{aa}));  % vitesse périgée [km/s]
        Va{aa} = sqrt(mu_e*(2/ra{aa}-1/a{aa}));  % vitesse apogée  [km/s]
        
        n{aa}  = sqrt(mu_e./a{aa}^3);        % Moyen mouvement   [rad/s]
        p{aa}  = a{aa}*(1-e{aa}^2);              % semilatus rectus  [km]
        
        T{aa}  = 2*pi/n{aa};                 % période          [s]
        h{aa}  = sqrt(p{aa}*mu_e);           % moment cinétique [km^2/s]
        
        h1{aa} = sin(i{aa})*sin(Ascension_droite_NA{aa});      % composant x du vecteur h
        h2{aa} = -sin(i{aa})*cos(Ascension_droite_NA{aa});     % composant y du vecteur h
        h3{aa} = cos(i{aa});                               % composant z du vecteur h
        
        n1{aa} = -h2{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant x de la ligne des noeuds
        n2{aa} =  h1{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant y de la ligne des noeuds
        n3{aa} = 0;                     % composant z de la ligne des noeuds
        
        N{aa}  = [n1{aa},n2{aa},n3{aa}];            % Vecteur ligne des noeuds
        
        
        
        
        fprintf('\n\n-------------- Paramètres orbitaux %d------------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        fprintf('Rayon du périgée numéro %d = %f Km \n',aa,rp{aa})
        fprintf('Altitude du périgée numéro %d = %f Km \n',aa,rp{aa}-R_e)
        fprintf('Rayon de l''apogée numéro %d = %f Km \n',aa,ra{aa})
        fprintf('Altitude de l''apogée numéro %d = %f Km \n\n',aa,ra{aa}-R_e)
        
        fprintf('Vitesse du périgée numéro %d = %f Km/s \n',aa,Vp{aa})
        fprintf('Vitesse à l''apogée numéro %d = %f Km/s \n\n',aa,Va{aa})
        
        fprintf('Période numéro %d = %f s \n',aa,T{aa})
        heures{aa}   = floor(T{aa}/3600);
        minutes{aa}  = floor((T{aa}-heures{aa}*3600)/60);
        secondes{aa} = floor(T{aa}-heures{aa}*3600-minutes{aa}*60);
        fprintf('Période numéro %d = %d h: %d m: %d s \n\n',aa,heures{aa},minutes{aa},secondes{aa});
        
        fprintf('Moyen mouvement numéro %d = %f rad/s \n',aa,n{aa})
        fprintf('Semilatus rectus numéro %d = %f Km \n\n',aa,p{aa})
        
        fprintf('Moment cinétique numéro %d = %f km^2/s \n',aa,h{aa})
        fprintf('Vecteur du moment (normalisé) numéro %d = [%f,%f,%f] \n',aa,h1{aa},h2{aa},h3{aa})
        fprintf('Vecteur ligne des noeuds numéro %d= [%f,%f,%f] \n',aa,n1{aa},n2{aa},n3{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n')
        
        
        
        fprintf('Rentrez le nombre d''orbite minimum souhaité :\n')
        norb{aa}=-1;
        while norb{aa}<=0 || floor(norb{aa})~=norb{aa}
            norb{aa} = input('');
        end
        
        
        delta_T{aa}=norb{aa}*T{aa};
        
        
        
        
        cos_E0{aa} = (e{aa}+cos(Nu{aa}))./(1+e{aa}.*cos(Nu{aa}));
        sin_E0{aa} = (sqrt(1-e{aa}^2).*sin(Nu{aa}))./(1+e{aa}.*cos(Nu{aa}));
        E0{aa} = atan2(sin_E0{aa},cos_E0{aa}); % Excentricité de l'anomalie (on parle de l'excentricité initiale) [rad]
        
        if (E0{aa}<0) 	% E0 appartient à [0,2pi]
            E0{aa}=E0{aa}+2*pi;
        end
        
        % pour le périgée qui se situr dans la suite :
        ellipse_per{aa}=(-E0{aa}+e{aa}.*sin(E0{aa}))./n{aa};
        
        
    elseif mmm{aa}==2
        %trajectoire parabolique
        
  
        fprintf('\nRentrez l''Ascension droite du nœud ascendant entre [0,360[ en degré:\n')
        Ascension_droite_NA{aa} = -1;
        while Ascension_droite_NA{aa}<0 || Ascension_droite_NA{aa}>360 
            Ascension_droite_NA{aa} = input('');
        end
        
        fprintf('Rentrez l''argument du périgée entre [0,360[ en degré:\n')
        Arg_p{aa} = -1;
        while Arg_p{aa}<0 || Arg_p{aa}>360
            Arg_p{aa}=input('');
        end
        
        fprintf('Rentrez l''anomalie vraie (pour t=0) entre [0,360[ en degré:\n')
        Nu{aa} = -1;
        while Nu{aa}<0 || Nu{aa}>360
            Nu{aa} = input('');
        end
        
        fprintf('Rentrez l''inclinaison entre [-90, 90] en degré:\n')
        i{aa} = -91;
        while i{aa}<-90 || i{aa}>90
            i{aa} = input('');
        end
        
        fprintf('Rentrez la distance au périgée de la terre (>%f) en Km:\n',R_e)
        rp{aa} = 0;
        while rp{aa}<R_e
            rp{aa} = input('');
        end
        
        
        
        fprintf('\n\n---------- Caracteristiques de l''orbite %d--------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        
        fprintf('Ascension droite du nœud ascendant numéro %d = %f deg \n',aa,Ascension_droite_NA{aa})
        fprintf('Argument au périgée numéro %d = %f deg \n',aa,Arg_p{aa})
        fprintf('Anomalie vraie au départ numéro %d = %f deg \n',aa,Nu{aa})
        fprintf('Inclinaison numéro %d = %f deg \n',aa,i{aa})
        fprintf('Distance au périgée numéro %d= %f en Km \n',aa,rp{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n')
        
        
        
        [Ascension_droite_NA_deg{aa},Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa}]=...
            radian_converter(Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa});
        
        
        
        
        %% Calcul des paramètres orbitaux
        
        
        Vp{aa} = sqrt(2*mu_e/rp{aa});           % vitesse périgée [km/s]
        
        p{aa}  = 2*rp{aa};                      % semilatus rectus  [km]
        h{aa}  = sqrt(2*mu_e*rp{aa});           % moment cinétique [km^2/s]
        
        h1{aa} = sin(i{aa})*sin(Ascension_droite_NA{aa});      % composant x du vecteur h
        h2{aa} = -sin(i{aa})*cos(Ascension_droite_NA{aa});     % composant y du vecteur h
        h3{aa} = cos(i{aa});                               % composant z du vecteur h
        
        n1{aa} = -h2{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant x de la ligne des noeuds
        n2{aa} =  h1{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant y de la ligne des noeuds
        n3{aa} = 0;                     % composant z de la ligne des noeuds
        
        N{aa}  = [n1{aa},n2{aa},n3{aa}];            % Vecteur ligne des noeuds
        
        e{aa}=1;
        
        fprintf('\n\n-------------- Paramètres orbitaux %d------------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        
        
        fprintf('Rayon du périgée numéro %d = %f Km \n',aa,rp{aa})
        fprintf('Altitude du périgée numéro %d = %f Km \n',aa,rp{aa}-R_e)
        
        fprintf('Vitesse du périgée numéro %d = %f Km/s \n',aa,Vp{aa})
        fprintf('Semilatus rectus numéro %d = %f Km \n\n',aa,p{aa})
        
        fprintf('Moment cinétique numéro %d = %f km^2/s \n',aa,h{aa})
        fprintf('Vecteur du moment (normalisé) numéro %d = [%f,%f,%f] \n',aa,h1{aa},h2{aa},h3{aa})
        fprintf('Vecteur ligne des noeuds numéro %d= [%f,%f,%f] \n',aa,n1{aa},n2{aa},n3{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n\n')
        
        
        
        %% Calcul de l'anomalie excentrique et de l'anomalie moyenne
        
        tp{aa}=input('Rentrez le temps au périgée : ');
        
        
    elseif mmm{aa}==3
        %trajectoire hyperbolique
        
        fprintf('\nRentrez l''Ascension droite du nœud ascendant entre [0,360[ en degré:\n')
        Ascension_droite_NA{aa} = -1;
        while Ascension_droite_NA{aa}<0 || Ascension_droite_NA{aa}>360 
            Ascension_droite_NA{aa} = input('');
        end
        
        fprintf('Rentrez l''argument du périgée entre [0,360[ en degré:\n')
        Arg_p{aa} = -1;
        while Arg_p{aa}<0 || Arg_p{aa}>360
            Arg_p{aa}=input('');
        end
        
        fprintf('Rentrez l''anomalie vraie (pour t=0) entre [0,360[ en degré:\n')
        Nu{aa} = -1;
        while Nu{aa}<0 || Nu{aa}>360
            Nu{aa} = input('');
        end
        
        fprintf('Rentrez l''inclinaison entre [-90, 90] en degré:\n')
        i{aa} = -91;
        while i{aa}<-90 || i{aa}>90
            i{aa} = input('');
        end
        
        fprintf('Rentrez le demi-grand axe de la trajectoire (>%f) en Km:\n',R_e)
        a{aa} = 0;
        while a{aa}<R_e
            a{aa} = input('');
        end
        
        
        ecc_max = 1+R_e/a{aa};
        fprintf('Rentrez l''excentricité. Il ne faut pas être en dessous de %f :\n', ecc_max)
        e{aa} = -1;
        while e{aa}<=ecc_max
            e{aa} = input('');
        end
        
        
        Nu_inf{aa}=acos(1/e{aa});
        
        
        
        fprintf('\n\n---------- Caracteristiques de l''orbite %d--------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        
        fprintf('Ascension droite du nœud ascendant numéro %d = %f deg \n',aa,Ascension_droite_NA{aa})
        fprintf('Argument au périgée numéro %d = %f deg \n',aa,Arg_p{aa})
        fprintf('Anomalie vraie au départ numéro %d = %f deg \n',aa,Nu{aa})
        fprintf('Inclinaison numéro %d = %f deg \n',aa,i{aa})
        fprintf('Anomalie vraie à l''infini numéro %d= %f deg \n',aa,Nu_inf{aa}*180/pi)
        fprintf('Demi-grand axe numéro %d = %f en Km \n',aa,a{aa})
        fprintf('Excentricité numéro %d = %f deg \n',aa,e{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n')
        
        
        
        
        [Ascension_droite_NA_deg{aa},Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa}]=...
            radian_converter(Ascension_droite_NA{aa},Arg_p{aa},Nu{aa},i{aa});
        
        
        
        
        
        
        %% Calcul des paramètres orbitaux
        
        
        
        rp{aa} = a{aa}*(e{aa}-1);                % rayon périgée [km]
        ra{aa} = -a{aa}*(e{aa}+1);                % rayon apogée  [km]
        
        
        Vp{aa} = sqrt(mu_e*(2/rp{aa}+1/a{aa})); % vitesse périgée [km/s]
        Va{aa} = sqrt(mu_e*(2/ra{aa}+1/a{aa})); % vitesse apogée  [km/s]
        V_inf{aa} = sqrt(mu_e/a{aa}); % vitesse a l'infini  [km/s]
        
        
        n{aa}  = sqrt(mu_e./a{aa}^3);         % Moyen mouvement   [rad/s]
        p{aa}  = a{aa}*(e{aa}^2-1);               % semilatus rectus  [km]
        
        
        h{aa}  = sqrt(p{aa}*mu_e);           % moment cinétique [km^2/s]
        
        h1{aa} = sin(i{aa})*sin(Ascension_droite_NA{aa});      % composant x du vecteur h
        h2{aa} = -sin(i{aa})*cos(Ascension_droite_NA{aa});     % composant y du vecteur h
        h3{aa} = cos(i{aa});                               % composant z du vecteur h
        
        n1{aa} = -h2{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant x de la ligne des noeuds
        n2{aa} =  h1{aa}/(sqrt(h1{aa}^2+h2{aa}^2)); % composant y de la ligne des noeuds
        n3{aa} = 0;                     % composant z de la ligne des noeuds
        
        N{aa}  = [n1{aa},n2{aa},n3{aa}];            % Vecteur ligne des noeuds
        
        
        
        
        fprintf('\n\n-------------- Paramètres orbitaux %d------------- \n',aa)
        fprintf('-------------------------------------------------- \n')
        fprintf('Rayon du périgée numéro %d = %f Km \n',aa,rp{aa})
        fprintf('Altitude du périgée numéro %d = %f Km \n',aa,rp{aa}-R_e)
        fprintf('Rayon de l''apogée numéro %d = %f Km \n',aa,ra{aa})
        fprintf('Altitude de l''apogée numéro %d = %f Km \n\n',aa,ra{aa}-R_e)
        
        fprintf('Vitesse du périgée numéro %d = %f Km/s \n',aa,Vp{aa})
        fprintf('Vitesse à l''apogée numéro %d = %f Km/s \n\n',aa,Va{aa})
        fprintf('Vitesse à l''infini numéro %d = %f Km/s \n\n',aa,V_inf{aa})
        
        fprintf('Moyen mouvement numéro %d = %f rad/s \n',aa,n{aa})
        fprintf('Semilatus rectus numéro %d = %f Km \n\n',aa,p{aa})
        
        fprintf('Moment cinétique numéro %d = %f km^2/s \n',aa,h{aa})
        fprintf('Vecteur du moment (normalisé) numéro %d = [%f,%f,%f] \n',aa,h1{aa},h2{aa},h3{aa})
        fprintf('Vecteur ligne des noeuds numéro %d= [%f,%f,%f] \n',aa,n1{aa},n2{aa},n3{aa})
        
        fprintf('-------------------------------------------------- \n')
        fprintf('-------------------------------------------------- \n\n')
        
    
    
        %% Calcul de l'anomalie excentrique et de l'anomalie moyenne
        
        
        sinF0{aa}=sqrt(e{aa}+1)+sqrt(e{aa}-1)*tan(Nu_inf{aa}/2);
        cosF0{aa}=sqrt(e{aa}+1)-sqrt(e{aa}-1)*tan(Nu_inf{aa}/2);
        
        F0{aa}=log(sinF0{aa}/cosF0{aa});
        
        if (F0{aa}<0) 	% E0 appartient à [0,2pi]
            F0{aa}=F0{aa}+2*pi;
        end
        
        tp{aa} = (e{aa}.*sinh(F0{aa})-F0{aa})./n{aa}; %+t0    % temps de passage périgée        [s]
    end
end


%% Angle Horaire de Greenwich

fprintf('\n\nAngle holraire du Méridien de Greenwich\n\n')
fprintf('Connaissez-vous l''Epoch ou le GHA ?\n');
fprintf('1- Epoch\n');
fprintf('2- GHA\n');
fprintf('Rentrez 1 ou 2\n')


zzz=0;
while zzz~=1 && zzz~=2
    zzz=input('');
end
    

if zzz==1
    GMST=GMST_converter();
    greenwich0=Ascension_droite_NA_deg{1} - GMST;
elseif zzz==2
    greenwich0 = input('\n\n Rentrez la longitude de Greenwich (le GHA) : ');
else
    fprintf('Il fallait choisir 1 ou 2')
    fprintf('Veuillez réessayer')
    return
end


%% Défintion du temps d'étude

delta_tmax=0;
tp_min=inf;
tp_max=0;
ellipse_per_max=0;

for aa=1:k
    if length(delta_T{aa})==1
        delta_tmax=max(delta_tmax,delta_T{aa});
    end
    
    if length(tp{aa})==1
        tp_max=max(tp_max,tp{aa});
        tp_min=min(tp_min,tp{aa});
    end
    
    if length(ellipse_per{aa})==1
        ellipse_per_max=max(ellipse_per_max,ellipse_per{aa});
    end
end

if tp_min==inf
    tp_min=0;
end

if tp_max==0
    tp_max=1;
end

delta_tmax=max(delta_tmax,ellipse_per_max);


fprintf('Vous allez devoir rentrer le temps initial et le temps final d''étude\n')
fprintf('Voici les contraites temporelles : \n')
fprintf('Il faut que le temps initial soit inférieure à %f\n',tp_min)
fprintf('Il faut que le temps final soit supérieure à %f\n',tp_max)
fprintf('Il faut que l''intervalle soit supérieure à %f\n\n',delta_tmax)

fprintf('Veuillez rentrer le temps initial\n')
t0=tp_min+1;
while t0>tp_min
    t0=input('');
end

    
fprintf('\nVeuillez rentrer le temps final\n')
tf=tp_max-1;
while (tf<tp_max) || (delta_tmax>tf-t0)
    tf=input('');
end


fprintf('La longeur de l''intervalle de temps est de %f s\n',tf-t0)


fprintf('\nVeuillez rentrer le pas\n')
pas=2*(tf-t0);
while pas>tf-t0
    pas=input('');
end

t=t0:pas:tf;



fprintf('\n\n------------------ Temps d''étude ----------------- \n')
fprintf('-------------------------------------------------- \n')
fprintf('Temps initial = %f s\n',t0)
fprintf('Temps final = %f s\n\n',tf)

fprintf('\nPas de discretisation = %f s\n',pas)
fprintf('Nombres de points = %f \n',length(t))

fprintf('-------------------------------------------------- \n')
fprintf('-------------------------------------------------- \n')



%% Définition rotation de la Terre
rot_earth  = Omega_e.*(t-t0)+greenwich0;  % GHA at the time t [rad]



%% Calcul de la position du satellite
M=cell(1,k);
E=cell(1,k);
sin_Nu=cell(1,k);
cos_Nu=cell(1,k);
v=cell(1,k);
theta=cell(1,k);
r=cell(1,k);

F=cell(1,k);
a_nu=cell(1,k);
b_nu=cell(1,k);


xp=cell(1,k);
yp=cell(1,k);

xs=cell(1,k);
ys=cell(1,k);
zs=cell(1,k);


rs=cell(1,k);
LatSSP=cell(1,k);
LongSSP=cell(1,k);

xSSP=cell(1,k);
ySSP=cell(1,k);
zSSP=cell(1,k);

for aa=1:k
    if mmm{aa}==1
        % ellipse
        
        
        %% Calcul de l'anomalie excentrique et de l'anomalie moyenne
        
        tp{aa} = (-E0{aa}+e{aa}.*sin(E0{aa}))./n{aa}+t0;    % temps de passage périgée        [s]
        M{aa}  = n{aa}.*(t-tp{aa});                 % anomalie moyenne                [rad]
        
        
        %% Mk = Ek - e*sin(Ek);
        
        E{aa} = zeros(size(t,2),1); % Vecteur excentricité anomalie
        for j=1:size(t,2)
            E{aa}(j) = anom_ecc(M{aa}(j),e{aa});  % excentricité anomalie  [rad]
        end
        
        
        %%  Anomalie vraie, Argument du périastre, Rayon
        
        sin_Nu{aa} = (sqrt(1-e{aa}.^2).*sin(E{aa}))./(1-e{aa}.*cos(E{aa}));
        cos_Nu{aa} = (cos(E{aa})-e{aa})./(1-e{aa}.*cos(E{aa}));
        v{aa} = atan2(sin_Nu{aa},cos_Nu{aa});         % Anomalie vraie (Nu en fonction du temps) [rad]
        
        theta{aa} = v{aa} + Arg_p{aa};                % Argument du périastre [rad]
        
        r{aa} = (a{aa}.*(1-e{aa}.^2))./(1+e{aa}.*cos(v{aa})); % Rayon  [km]
        
        
    elseif mmm{aa}==2
        % parabole
        
        
        
        %% Calcul de l'anomalie moyenne
        
        M{aa}  = (mu_e^2)/(h{aa}^3)*(t-tp{aa});               % anomalie moyenne                [rad]
        
        
        %% Equation de BARKER
        
        v{aa} = zeros(size(t,2),1);
        for j=1:size(t,2)
            A = nthroot(3*M{aa}(j)-sqrt(9*(M{aa}(j)^2)+1),3);
            B = nthroot(3*M{aa}(j)+sqrt(9*(M{aa}(j)^2)+1),3);
            
            v{aa}(j) = 2*atan(A+B);  % anomalie vraie  [rad]
        end
        
        
        %%  Argument du périastre, Rayon
        
        theta{aa} = v{aa} + Arg_p{aa};                % Argument du périastre [rad]
        
        r{aa} = ((h{aa}^2)/mu_e)./(1+cos(theta{aa}-Arg_p{aa})); % Rayon  [km]
        
        
        
        
        
    elseif mmm{aa}==3
        % hyperbole
        
        
        M{aa}  = n{aa}.*(t-tp{aa});                % anomalie moyenne                [rad]
        
        
        %% Mk = e*sinh(Fk) - Fk ;
        
        
        F{aa} = zeros(size(t,2),1); % Vecteur excentricité anomalie
        for j=1:size(t,2)
            F{aa}(j) = anom_hyperbolic(M{aa}(j),e{aa});  % excentricité anomalie  [rad]
        end
        
        %%  Anomalie vraie, Argument du périastre, Rayon
        
        a_nu{aa} = sqrt((e{aa}+1)/(e{aa}-1));
        b_nu{aa} =tanh(F{aa}/2);
        v{aa} = 2*atan(a_nu{aa}*b_nu{aa});         % Anomalie vraie (Nu en fonction du temps) [rad]
        
        theta{aa} = v{aa} + Arg_p{aa};                % Argument du périastre [rad]
        
        r{aa} = (a{aa}.*((e{aa}.^2)-1))./(1+e{aa}.*cos(v{aa})); % Rayon  [km]
        
        
    end
    
    
    
    
    %% Coordonnées du Satellite
    % "Inertial" reference system ECI (Earth Centered Inertial)
    
    xp{aa} = r{aa}.*cos(theta{aa});      % x position (direction noeud)   [km]
    yp{aa} = r{aa}.*sin(theta{aa});      % y position (perpendiculaire x) [km]
    
    % Coordonnées (x,y,z) dans le repère ECI [Km]
    xs{aa} = xp{aa}.*cos(Ascension_droite_NA{aa})-yp{aa}.*cos(i{aa}).*sin(Ascension_droite_NA{aa});
    ys{aa} = xp{aa}.*sin(Ascension_droite_NA{aa})+yp{aa}.*cos(i{aa}).*cos(Ascension_droite_NA{aa});
    zs{aa} = yp{aa}.*sin(i{aa});
    
    rs{aa} = p{aa}./(1+e{aa}.*cos(theta{aa}-Arg_p{aa}));       % norm radius sat  [km]
    
    
    %% Greenwich hour angle (GHA)
    
    for j=1:size(t,2)
        if rot_earth(j) < (-pi)
            nn = ceil(rot_earth(j)/(-2*pi));
            rot_earth(j) = rot_earth(j) + nn*2*pi;
        elseif rot_earth(j) > (pi)
            nn = fix(rot_earth(j)/(2*pi));
            rot_earth(j) = rot_earth(j) - nn*2*pi;
        end
    end
    
    LatSSP{aa}    = asin(sin(i{aa}).*sin(theta{aa}));          % Latitude Satellite   [rad]
    LongSSP{aa}    = atan2(ys{aa}./rs{aa},xs{aa}./rs{aa})-rot_earth';   % Longitude Satellite  [rad]
    
    xSSP{aa} = R_e.*cos(LatSSP{aa}).*cos(LongSSP{aa})+1;         % x Satellite        [km]
    ySSP{aa} = R_e.*cos(LatSSP{aa}).*sin(LongSSP{aa})+1;         % y Satellite        [km]
    zSSP{aa} = R_e.*sin(LatSSP{aa})+0.1;                     % z Satellite        [km]
    
end



%% Vues pour le plot


str = 0;
fprintf(' ');
fprintf(' \n\n********** Options vue 3D  *********** \n');
fprintf(' ----------- Point de vue ------------ \n');
fprintf(' 1) Moment cinétique \n');
fprintf(' 2) Equateur \n');
fprintf(' 3) Pole Nord \n');
fprintf(' 4) Suivre le satellite dans toutes les dimensions (semble sans mouvement) \n');
fprintf(' 5) Suit le satellite dans la Longitude \n');
fprintf(' 6) Suit le méridien de Greenwich \n');
fprintf('    autre \n');
pw  = input(' Rentrez 1, 2, 3, 4, 5, 6 ou autre : ','s');
if (exist('pw','var')==0)
    pw = '0';
end

fprintf(' \n\n---------------- Zoom ---------------- \n');
fprintf(' 1) Aucun\n');
fprintf(' 2) Zoom Normal \n');
fprintf(' 3) Zoom Maximal \n');
while ((str~=1)&&(str~=2)&&(str~=3))
    str = input(' Rentrez 1, 2 ou 3 : ');
end




%% FIGURE

% Charger image Terre
cdata = imread(image_file);

a_max=0;
rp_max=0;
for aa=1:k
    if length(a{aa})==1
        a_max=max(a_max,a{aa});
    end
    
    if length(rp{aa})==1
        rp_max=max(rp_max,rp{aa});
    end
end

aa_max=0;
for aa=1:k
    if a_max==a{aa}
        aa_max=aa;
        break
    end
end

for aa=1:k
    if length(a_max)==0  %#ok<ISMT>
        a_max=rp_max;
        
    elseif length(a_max)==1 && length(rp_max)==1
        a_max=max(a_max,rp_max);
    end
    
end


%% VIEW 3D ORBIT

X=cell(1,k);
XS=cell(1,k);

% Affichage equateur
angle_eq = linspace(0,2*pi,361);
xeq      = (R_e*1.0001).*cos(angle_eq);
yeq      = (R_e*1.0001).*sin(angle_eq);
zeq      = zeros(1,size(angle_eq,2));

% Affichage
figure('Color','k','units','normalized','outerposition',[0 0 1 1])
% Definition of axes necessary to rotate the object to whom are reported
ax = axes('XLim',[-a_max a_max],'YLim',[-a_max a_max],'ZLim',[-a_max a_max],'Color','k');

% Centre de la Terre
[x,y,z] = ellipsoid(0, 0, 0, R_e, R_e, R_e, 30);
hold on;
switch lower(pw)
    case '1'
        view([h1{aa_max} h2{aa_max} h3{aa_max}]);
    case '2'
        view(0,0);
    case '3'
        view(120,90);
end

% plot the surface ellipsoid
globe = surface(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
% applying the properties of the surface setting up the image_file
set(globe, 'FaceColor', 'texturemap', 'CData', cdata,'FaceAlpha', 0.9, 'EdgeColor', 'none');

grid on;
axis equal;

% empiric formulation to apply the adequate zoom
if ((a_max/R_e<13.3)&&(str==2))
    if (a_max/R_e>9)
        zoom(13.3-a_max/R_e);
    elseif (a_max/R_e>6.8)
        zoom(10.7-a_max/R_e);
    elseif (a_max/R_e>5)
        zoom(9.6-a_max/R_e);
    elseif (a_max/R_e>3.5)
        zoom(6.1-a_max/R_e);
    elseif (a_max/R_e>2.7)
        zoom(5.4-a_max/R_e);
    elseif (a_max/R_e>2.2)
        zoom(5.0-a_max/R_e);
    elseif (a_max/R_e>1.5)
        zoom(3.3-a_max/R_e);
    else
        zoom(2.7-a_max/R_e);
    end
elseif (str==3)
    zoom(a_max^1.03/R_e);
elseif (str==1)
    zoom(1);
end

% limit of axes
axis([-a_max a_max -a_max a_max -a_max a_max]);
% Declaration of variable that indicates a ratio of connection with axes
hg = hgtransform('Parent',ax);
% Set the globe to depend on axes (if axes rotate thet also globe rotates)
set(globe,'Parent',hg);
set(gcf,'Renderer','opengl');

% Inertial reference system, equator and nodes' line
plot3([0,2*a_max],[0,0],[0,0],'-.w','LineWidth',1);    % X = Aries' direction
plot3([0,0],[0,2*a_max],[0,0],'-.w','LineWidth',1);    % Y
plot3([0,0],[0,0],[0,2*a_max],'-.w','LineWidth',1);    % Z
text(2*a_max+120,10,0,texlabel('gamma'),'Color','w','FontSize',18);
text(10,2*a_max+120,0,texlabel('Y'),'Color','w','FontSize',10);
text(0,0,2*a_max+140,texlabel('Z'),'Color','w','FontSize',10);
plot3(xeq,yeq,zeq,'--w','LineWidth',1);                           % Equator
plot3([0,2*a_max*n1{aa_max}],[0,2*a_max*n2{aa_max}],[0,n3{aa_max}],'-.y','LineWidth',1.5);     % Nodes' Line
text(2*a_max*n1{aa_max}-140,2*a_max*n2{aa_max}+140,0,texlabel('RAAN'),'Color','y','FontSize',8);

u = 0;                                   % initialization of counter
for kk=2:size(t,2)                        % cycle for setting the view required
    if strcmp(pw,'4')
        view([xs{1}(kk),ys{1}(kk),zs{1}(kk)]);
    elseif strcmp(pw,'5')
        view([xs{1}(kk),ys{1}(kk),2]);
    elseif strcmp(pw,'6')
        view(rot_earth(kk)*180/pi+89,2);
    elseif strcmp(pw,'0')
        view(120,45);
    end                       % Property of rotation for axes
    
    
    Rz = makehgtform('zrotate',rot_earth(kk));
    
    for aa=1:k
        
        X{aa} = plot3(xSSP{aa}(kk),ySSP{aa}(kk),zSSP{aa}(kk),'--go','LineWidth',2,...         % Plot Sub-Satellite-Point
            'MarkerSize',2.5,'MarkerEdgeColor',colors{aa},'MarkerFaceColor',colors{aa});
        
        set(X{aa},'Parent',hg);           % Set the SSP to depend on axes (if axes rotate thet also SSP rotates)
        set(hg,'Matrix',Rz);          % Application of the rotation of angle wE*time
        drawnow                       % draw the rotation
        
        
        u = u+1;
        XS{aa}(kk) = plot3(xs{aa}(kk),ys{aa}(kk),zs{aa}(kk),'rp','LineWidth',0.9,...    % Plot the sudden position of the satellite
            'MarkerSize',10,'MarkerEdgeColor',colors{aa},'MarkerFaceColor',colors{aa});
        pause(0.0001);
        
        if (u>1)
            delete(XS{aa}(kk-1));               % delete last position of the satellite and replace it with ...
            XS{aa}(kk-1) = plot3(xs{aa}(kk-1),ys{aa}(kk-1),zs{aa}(kk-1),'--x','MarkerEdgeColor',colors{aa},'LineWidth',0.8,'MarkerSize',2); % ... symbol of the orbit
        end
    end
end




for j=1:size(t,2)
    for aa=1:k
        if LongSSP{aa}(j) < (-pi)             % Longitude entre [-180,180] deg
            %nn = fix(-LongSSP(j)/(pi));
            nn = 1;
            LongSSP{aa}(j) = LongSSP{aa}(j) + 2*nn*pi;
        elseif LongSSP{aa}(j) > (pi)
            %nn = fix(LongSSP(j)/(pi));
            nn = 1;
            LongSSP{aa}(j) = LongSSP{aa}(j) - 2*nn*pi;
        end
    end
end

for aa=1:k
    LongSSP{aa} = LongSSP{aa}.*180/pi;             % Longitude of SSP [deg]
    LatSSP{aa}  = LatSSP{aa}.*180/pi;              % Latitude  of SSP [deg]
end

figure('Name','Planisphere','units','normalized','outerposition',[0 0 1 1])
hold on;
set(gca,'XTick',[-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
    15 30 45 60 75 90 105 120 135 150 165 180],'XTickMode','manual');
set(gca,'YTick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'YTickMode','manual');

imagesc([-180,180],[90,-90],cdata);

grid on;
axis([-180 180 -90 90]);
plot([-180 180],[0 0],'--k','LineWidth',1.5);  % Equator


for j=2:size(t,2) % cycle for plotting all orbits
    for aa=1:k
        %aa
        plot(LongSSP{aa}(j),LatSSP{aa}(j),'Color',colors{aa},'MarkerSize',1,'MarkerFaceColor',[0.8 0.2 0.1],...
            'MarkerEdgeColor',colors{aa}); % a modifier
        pause(0.05);
        if (j>1)
            if (abs(LongSSP{aa}(j)-LongSSP{aa}(j-1))<80)
                plot([LongSSP{aa}(j-1),LongSSP{aa}(j)],[LatSSP{aa}(j-1),LatSSP{aa}(j)],'-','Color',colors{aa},...
                    'LineWidth',1.5); % a modifier
            end
        end
    end
end