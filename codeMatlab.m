clear;close all;

%Chargement des données de Vent :
load("data/vent_global.mat");

%Chargmement des données géographiques
num = readtable("data/Carte_globale_tous_pays_modif.xlsx");
pluieTemp = table2array(num(:,21));
pays = table2array(num(:,14));

%Caractéristiques des polluants et de leurs sites pollueurs
Polluants = ["PM106","PM109","PM1014","PM1015"];
CDUh = ones(4,1)*10;
Solubilite = zeros(4,1);
Tsites = [19,33;54,17;37,43;61,37];
Polsites = [24 0 0 0 ; 0 9.176027546 0 0 ; 0 0 90.93535659/2 0 ; 0 0 0 9.176027546];

%Nombre de site
dim = size(Tsites);
Nsites = dim(1);

%Dimensions spatiales et temporelles domaine réel
deltaT = 3*60;
deltaTnorm = 3*60*60;
deltaX = 35000;
deltaXnorm = 70000;

%Variable de la loi de Fick
rho = 1.293;
D = 500;

%Dimensions spatiales et temporelles domaine numérique
dim = size(uwtot);
Nt = dim(3);
Nx = dim(2)*deltaXnorm/deltaX;
NxVent = dim(2);
Ny = dim(1)*deltaXnorm/deltaX;
NyVent = dim(1);
dim = size(Polluants);
Np = dim(2);
clear dim;

%Initialisation des matrices
C = zeros(Nx+2,Ny+2,Np,2);%Quantité de polluant dans l'air
totalite = zeros(Nx+2,Ny+2,Np);%Moyenne de présence de polluant dans l'air sur la totalité de la durée de simulation
Msol = zeros(Nx+2,Ny+2,Np);%Quantité de polluant au sol
pluiev2 = zeros(NxVent+2,NyVent+2,Np);% Quantité de pluie (mm/an)
paysR = zeros(NxVent+2,NyVent+2,Np);%Pays

%Définition des pays et de la pluie
for x=1:NxVent
    for y=1:NyVent
        n = findEmpl(num,x,y);
        if(n < 0)
            pluiev2(x+1,y+1) = 0;
            paysR(x+1,y+1) = "Neutral";
        else
            pluiev2(x+1,y+1)=pluieTemp(n);
            paysR(x+1,y+1) = pays(n);
        end
    end
end

prc = 0;

%Simulation
for t=2:(Nt*(deltaTnorm/deltaT))
    for xt=1:(Nx+2)
        for yt=1:(Ny+2)
            for i=1:Np

                if(xt > 1 && yt > 1 && xt <= Nx+1 && yt <= Ny+1)
                    
                    %Définition des positions x,y dans la matrice des vents
                    x =round((xt-1)*deltaX/deltaXnorm)+1;
                    y =round((yt-1)*deltaX/deltaXnorm)+1;
                    if x > NxVent
                        x = NxVent;
                    end
                    if y > NyVent
                        y = NyVent;
                    end

                    %Source de pollution
                    S = 0;
                    for sites=1:Nsites
                        if x == Tsites(sites,1) && y == Tsites(sites,2)
                            S = Polsites(sites,i)*deltaT/(365*24*60*60);
                        end
                    end
                    
                    %Définition de la valeur du vent
                    tcorr = 1+round( t*(deltaT/deltaTnorm) );
                    if(tcorr > 120)
                        tcorr=120;
                    end
                    uxu = uwtot(y,x,tcorr);
                    uyu = vwtot(y,x,tcorr);
                     
                    
                    %Equation différentielle discrétisée
                    C(xt,yt,i,2) = C(xt,yt,i,1)+(deltaT/rho)*( D/(deltaX*deltaX)*(C(xt+1,yt,i,1)+C(xt-1,yt,i,1)+C(xt,yt+1,i,1)+C(xt,yt-1,i,1)-4*C(xt,yt,i,1)) + S - uxu*(C(xt+1,yt,i,1)-C(xt-1,yt,i,1))/(2*deltaX) - uyu*(C(xt,yt+1,i,1)-C(xt,yt-1,i,1))/(2*deltaX) );
                    if C(xt,yt,i,2) < 0
                        C(xt,yt,i,2) = 0;
                    end
                    
                    %Dépot de polluant au sol
                    MmaxDissouePluie = exp(-500/pluiev2(x,y))*C(xt,yt,i,2)*0.01;
                    MdeposerPluie = min(MmaxDissouePluie,C(xt,yt,i,2));
                    Msol(xt,yt,i) = Msol(xt,yt,i) + 0.1*MdeposerPluie;
                    C(xt,yt,i,2) = C(xt,yt,i,2) - 0.1*MdeposerPluie;

                    
                else
                    %Conditions limites
                    if(xt == 1)
                        C(xt,yt,i,2) = C(xt+1,yt,i,1);
                    elseif(xt == Nx+2)
                        C(xt,yt,i,2) = C(xt-1,yt,i,1);
                    elseif(yt == 1)
                        C(xt,yt,i,2) = C(xt,yt+1,i,1);
                    elseif(yt == Ny +2)
                        C(xt,yt,i,2) = C(xt,yt-1,i,1);
                    end
                    
                end

            end
        end
    end
    
    %Mise à jour des valeurs et moyenne
    C(:,:,:,1) = C(:,:,:,2);
    totalite(:,:,:) = totalite(:,:,:) + C(:,:,:,1)/(Nt*(deltaTnorm/deltaT));
    
    %Affichage 
    if mod(t,10)==0 || t == 2
        axis([0 Nx 0 Nx]);
        colorbar
        h = pcolor(transpose( C(:,:,1,2)+C(:,:,2,2)+C(:,:,3,2)+C(:,:,4,2)  ));
        set(h, 'EdgeColor', 'none');
        colorbar
        axis([0 Nx 0 Nx]);
        
        nom = string(t);
        while(strlength(nom) < 6)
            nom = "0"+nom;
        end
        
        %Pour sauvegarder les tracés :
        %saveas(h,nom+".png");
        
        title(num2str(round(t*(deltaT/deltaTnorm)*3))+"h");
    end
    pause(0.3);

    %Affichage du pourcentage d'avancement de la simulation :
    prcp = round(t*100/(Nt*(deltaTnorm/deltaT)));
    if prcp > prc
        prc = prcp
    end
    
end

vars = {"prc","D","h","i","MdeposerPluie","MmaxDissouePluie","n","nom","pays","pluieTemp","prcp","readme","S","starttime","t","tcorr","uxu","uyu","x","xt","y","yt"};
clear(vars{:});
clear vars;