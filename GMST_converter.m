function angle=GMST_converter(varargin)

%% Initialisation de la date et l'heure

jo = input('\n\nRentrez le numéro du jour entre 1 et 31 : ');
Day={'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',...
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',...
    '25', '26', '27', '28', '29', '30', '31'};
jour=char(Day(jo));


mo = input('Rentrez le numéro du mois entre 1 et 12 : ');
Month = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
mois = char(Month(mo));

annee=input('Rentrez l''année supérieure à 2000 : ','s');


ho = input('Rentrez l''heure entre 1 et 24 : ');
Hour={'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',...
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'};
heure = char(Hour(ho));

mi = input('Rentrez les minutes entre 0 et 60 : ');
Minutes={'00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',...
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',...
    '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36',...
    '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48',...
    '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60'};
minute = char(Minutes(mi+1));


so = input('Rentrez les secondes entre 0 et 60 : ');
seconde = char(Minutes(so+1));

fprintf('\n\n-------------- Affichage de l''Epoch -------------- \n')
fprintf('-------------------------------------------------- \n')
fprintf('Année : %s\n',annee)
fprintf('Mois : %s\n',mois)
fprintf('Jour : %s \n\n',jour)

fprintf('Heure : %s\n',heure)
fprintf('Minute : %s\n',minute)
fprintf('Seconde : %s\n',seconde)
fprintf('-------------------------------------------------- \n')
fprintf('-------------------------------------------------- \n')
%{
jour = '26';
mois = 'Dec';
annee = '2009';
heure = '12';
minute = '00';
seconde = '00';
%}

f = strcat(jour,'-',mois,'-',annee);
g=strcat(heure,':',minute,':',seconde);

str=[f ' ' g];

%% Calcul de l'angle

td = datetime(str);
q=juliandate(td);

angle=JD2GMST(q);

fprintf('\n\n-------------------- Résultats ------------------- \n')
fprintf('-------------------------------------------------- \n')
fprintf('date et heure : %s\n',td)
fprintf('Jour Julien : %f\n\n',q)
fprintf('Angle de l''Epoch : %f\n',angle)
fprintf('-------------------------------------------------- \n')
fprintf('-------------------------------------------------- \n')

end