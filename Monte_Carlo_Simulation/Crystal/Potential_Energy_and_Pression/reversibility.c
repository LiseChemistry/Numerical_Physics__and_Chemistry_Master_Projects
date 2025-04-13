#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

//-----------------------------{Définition des constantes}-----------------------------
const int nbx_particules = 100;   //Nombre de particules du système
#define Rc 2.5* sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  40               // Taille de la boîte
const int sigma = 1;                // sigma du lennard jones adimentionné
const double Lmoitie = L * 0.5;    // calcul de la valeur de la moitié de L
const double rho = nbx_particules / (L * L); // rho la densité volumique
const double epsilon = 1;          //le E0 du lennard jones adimentionné
const double V=L*L; 
const int dimension = 2; 

/*Calcul pour Rc pour eviter de les refaire, vu qu'ils sont constants*/
const double Rc2 = Rc * Rc;        
const double inv_Rc_2 = 1 / Rc2;
const double inv_Rc_4 = inv_Rc_2 * inv_Rc_2;
const double inv_Rc_6 = inv_Rc_2 * inv_Rc_2 * inv_Rc_2;
const double inv_Rc_10 = inv_Rc_6 * inv_Rc_4;
const double inv_Rc_12 = inv_Rc_6 * inv_Rc_6;

/*Calcul des terme à ajouter pour les grandeurs calculées, qui sont constants*/
const double terme_commun_frac = 2/5;
const double U_tail = pi * epsilon * nbx_particules * rho * (terme_commun_frac * inv_Rc_10 - inv_Rc_4);
const double U_decal  = 4 * (inv_Rc_12 - inv_Rc_6);
const double P_tail = 6*pi*epsilon*rho*rho*(terme_commun_frac*inv_Rc_10 - 0.5*inv_Rc_4);

//Je définis le delta max pour le placement aléatoire des particules
const double delta_max = 0.25*sigma; 

//-----------------------------{Définition de la structure pour la particule}-----------------------------
/**
 * Définition de la classe Particule avec les attribus.
 *
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule
 * @param vx vitesse selon x de la particule.
 * @param vy vitesse selon y de la particule.
 * @param fx force appliquée selon x sur la particule.
 * @param fy force appliquée selon y sur la particule.
 * @param actif  indique si la particule existe(sert pour les défauts)
 */
typedef struct {
    double x; 
    double y;
    double fx; 
    double fy; 
    int actif; //1=active et 0= inactive, permet d'inclure des défauts
} Particule;

//-----------------------------{Constructeur pour la particule}---------------------------------------
/**
 * @brief Initialise les attributs d'une instance de Particule.
 *
 * Cette fonction configure les coordonnées et les vitesses d'une particule 
 * en initialisant ses attributs avec les valeurs fournies.
 *
 * @param p Pointeur vers l'instance de Particule à initialiser.
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule.
 */
void init_Particule(Particule *p, double x, double y) {
    // Initialisation des coordonnées de la particule
    p->x = x;
    p->y = y;
}

void initialiser_Configuration_Cristalline(Particule tab_par[]){

    //Calculer le nombre de particules par ligne et par colonne
    int particules_par_ligne = (int) sqrt(nbx_particules); 
    int particules_par_colonne = (nbx_particules+particules_par_ligne-1)/particules_par_ligne; 

    //Espacement entre les particules
    double espacement_x = 4.0; 
    double espacement_y = 4.0; 

    //Calcul de la largeur et hauteur du réseau 
    double largeur_reseau = espacement_x*(particules_par_ligne-1); 
    double hauteur_reseau = espacement_y*(particules_par_colonne-1); 

    //Décalage pour centrer le réseau atour de (0,0) 
    double decalage_x = -largeur_reseau/2.0; 
    double decalage_y = -hauteur_reseau/2.0; 

    for (int i=0; i<= particules_par_ligne; i++){
        for (int j=0; j< particules_par_colonne; j++){

            int index = i*particules_par_colonne+j; 
            if (index<nbx_particules){

                // Position initiale centrée autour de (0, 0)
                tab_par[index].x = i * espacement_x + decalage_x;
                tab_par[index].y = j * espacement_y + decalage_y;

                tab_par[index].actif = 1;
            }
        }
    }
}



//-----------------------------{changement de configuration aleatoire}------------------------------------------
/** 
 * 
 * 
*/
int change_nouvelle_configuration(Particule tab_par[], int nbx_particules)
{
    // Générer un nombre entier aléatoire entre 0 et N exclus
    int index_particule = rand() % (nbx_particules);

    // Générer un rayon aléatoire entre 0 et delta_max
    double rayon = ((double)rand() / RAND_MAX) * delta_max;

    // Générer un angle aléatoire entre 0 et 2*pi (en radians)
    double angle = ((double)rand() / RAND_MAX) * 2 * M_PI;

    // Calculer les déplacements en x et y en utilisant le rayon et l'angle
    double delta_x = rayon * cos(angle);
    double delta_y = rayon * sin(angle);

    // Ajouter ces déplacements à la position actuelle de la particule
    tab_par[index_particule].x += delta_x;
    tab_par[index_particule].y += delta_y;

    //Application des conditions périodiques sur x : 
    if (tab_par[index_particule].x >= Lmoitie) {tab_par[index_particule].x -=L;}
    if (tab_par[index_particule].x < -Lmoitie) {tab_par[index_particule].x += L;}

    //Application des conditions périodiques sur y : 
    if (tab_par[index_particule].y >= Lmoitie) {tab_par[index_particule].y -=L;}
    if (tab_par[index_particule].y < -Lmoitie) {tab_par[index_particule].y += L;}

    return index_particule;
}

//-----------------------------{Calcul de Ep}-----------------------------------------

double calcul_energie_potentielle_terme_commun(Particule tab_par[], int index_modifie) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double U = 0.0; // energie potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L * L; //volume
    int dimension = 2;

    for (int i = 0; i < nbx_particules; i++) 
    {

        if (i!=index_modifie){

            /*condition qui verifie que le calcul se fait que par rapports aux partiucules et ignore les défauts*/
            if (tab_par[i].actif) 
            {
                for (int j = i + 1; j < nbx_particules; j++) 
                {

                    if (j!=index_modifie){

                        if (!tab_par[j].actif) continue; 


                    double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                    double y_sous = tab_par[j].y - tab_par[i].y; //

                    double abs_distx = fabs(x_sous);
                    double abs_disty = fabs(y_sous);

                    if (abs_distx > Lmoitie) {
                        x_sous = - (x_sous / abs_distx) * (L - abs_distx);
                    }

                    if (abs_disty > Lmoitie) {
                        y_sous = - (y_sous / abs_disty) * (L - abs_disty);
                    }
                        
                    /*Calcul de la distance entre deux particules*/
                    rij2 = x_sous * x_sous + y_sous * y_sous;   

                    /*Calcul sur les puissances de rij pour eviter l'appel de pow*/
                    inv_r2 = 1.0 / rij2;                  //   1/rij^2
                    inv_r2x3 = inv_r2 * inv_r2 * inv_r2; //    1/rji^6
                    inv_r2x6 = inv_r2x3 * inv_r2x3;     //    1/rij^12

                    /*Conditions qui permet de tronquer le potentiel à Rc*/
                    if (rij2 < Rc2) 
                    {
                        U += 4 * (inv_r2x6 - inv_r2x3) - U_decal; //calcul du potentiel entre deux particules
                    }
                    }
                }    
            }
        }
        
    }
    return U;    //calcul et mise à jour l'energie potentielle sans U_tail
}

double calcul_energie_potentielle_terme_modifie(Particule tab_par[], int index_modifie) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double U_index = 0.0; // energie potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L * L; //volume
    int dimension = 2;

    int i = index_modifie; 

    if (tab_par[i].actif){

        for (int j=0; j<nbx_particules; j++){
            if (!tab_par[j].actif) continue; 

            if(j!=i){
                
                double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                double y_sous = tab_par[j].y - tab_par[i].y; //

                double abs_distx = fabs(x_sous);
                double abs_disty = fabs(y_sous);

                if (abs_distx > Lmoitie) {
                    x_sous = - (x_sous / abs_distx) * (L - abs_distx);
                }

                if (abs_disty > Lmoitie) {
                    y_sous = - (y_sous / abs_disty) * (L - abs_disty);
                }
                
                /*Calcul de la distance entre deux particules*/
                rij2 = x_sous * x_sous + y_sous * y_sous;   

                /*Calcul sur les puissances de rij pour eviter l'appel de pow*/
                inv_r2 = 1.0 / rij2;                  //   1/rij^2
                inv_r2x3 = inv_r2 * inv_r2 * inv_r2; //    1/rji^6
                inv_r2x6 = inv_r2x3 * inv_r2x3;     //    1/rij^12

                /*Conditions qui permet de tronquer le potentiel à Rc*/
                if (rij2 < Rc2) 
                {
                    U_index += 4 * (inv_r2x6 - inv_r2x3) - U_decal; //calcul du potentiel entre deux particules
                }
            }
        }
    }

    return U_index;
}

//-----------------------------{calcul de boltzmann}------------------------------------------

double calcul_boltzmann(Particule tab_ini[], Particule tab_temp[], int index_modifie, double temperature)
{
    double proba = 0.0;

    double terme_commun = calcul_energie_potentielle_terme_commun(tab_ini, index_modifie) + U_tail;

    double energie_tab_ini = terme_commun + calcul_energie_potentielle_terme_modifie(tab_ini,index_modifie); 
    double energie_tab_temp = terme_commun + calcul_energie_potentielle_terme_modifie(tab_temp,index_modifie); 

    if (energie_tab_temp <= energie_tab_ini){
        proba = 1.0; 
    } else {
        proba = exp(-(energie_tab_temp-energie_tab_ini)/(temperature));
    }

    return proba; 
}



//-----------------------------{Enregistrer les positions}-------------------------------------------

void enregistrer_positions(FILE *fp, Particule *tab_par, int tentative, int nbx_particules){
    fprintf(fp, "%d\n", nbx_particules);
    fprintf(fp, "\n");
    for (int i=0; i<nbx_particules; i++){
        fprintf(fp, "A %f %f\n", tab_par[i].x, tab_par[i].y);
    }
}

void enregistrer_energie(FILE *fp, double energie, int cycle){
    fprintf(fp, "%d : %f\n", cycle, energie); 
}

void enregistrer_pression(FILE *fp, double pression, int cycle){
    fprintf(fp, "%d : %f\n", cycle,pression);
}

//-----------------------------{Vérifier les différences}-------------------------------------------

/**
 * Fonction qui permet de vérifier quelle particule a bougée entre deux configurations
 */
void verifier_differences(Particule tab_ini[], Particule tab_temp[], int nbx_particules) {
    for (int i = 0; i < nbx_particules; i++) {
        if (tab_ini[i].x != tab_temp[i].x || tab_ini[i].y != tab_temp[i].y) {
            printf("Particule %d a changé : \n", i);
            printf("  Initial : x = %f, y = %f\n", tab_ini[i].x, tab_ini[i].y);
            printf("  Temporaire : x = %f, y = %f\n", tab_temp[i].x, tab_temp[i].y);
        }
    }
}

//-----------------------------{Calcul_forces}--------------------------------------------------------
/**
 * @brief Calcule les forces entre toutes les paires de particules dans le système.
 *
 * Cette fonction parcourt le tableau des particules et calcule les forces d'interaction
 * entre toutes les paires, mettant ainsi à jour les composantes de force de chaque particule.
 * Le calcul utilise des conditions périodiques aux bords pour gérer les interactions
 * et contribue également au calcul du viriel, utilisé dans le calcul de la pression.
 *
 * @param tab_par Tableau contenant les particules du système, chacune avec ses
 *                attributs de position, de vitesse, et de force.
 * @param viriel Pointeur vers une variable où le viriel total du système sera stocké.
 *               Le viriel est une somme pondérée des forces entre particules et sert
 *               dans les calculs thermodynamiques (notamment la pression).
 *
 * @note La fonction suppose que `tab_par` est initialisé avec les particules en positions
 *       valides dans l'espace de simulation.
 */
double calcul_viriel(Particule tab_par[]) {
    // Initialiser le viriel (sert pour le calcul de la pression)
    double viriel = 0.0;

    // Réinitialiser les forces et l'énergie potentielle
    for (int i = 0; i < nbx_particules; i++) {
        tab_par[i].fx = 0.0;
        tab_par[i].fy = 0.0;
    }

    // Calcul des forces d'interaction entre les particules
    for (int i = 0; i < nbx_particules; i++) {
        if (!tab_par[i].actif) continue;

        for (int j = i + 1; j < nbx_particules; j++) {
            if (!tab_par[j].actif) continue;

            double x_sous = tab_par[i].x - tab_par[j].x;
            double y_sous = tab_par[i].y - tab_par[j].y;

            double abs_distx = fabs(x_sous);
            double abs_disty = fabs(y_sous);

            if (abs_distx > Lmoitie) {
                x_sous = - (x_sous / abs_distx) * (L - abs_distx);
            }

            if (abs_disty > Lmoitie) {
                y_sous = - (y_sous / abs_disty) * (L - abs_disty);
            }

            double rij2 = x_sous * x_sous + y_sous * y_sous;

            if (rij2 < Rc2) {
                double r6 = rij2 * rij2 * rij2;
                double r12 = r6 * r6;
                double inv_r2 = 1.0 / rij2;
                double inv_r6 = 1.0 / r6;
                double inv_r12 = 1.0 / r12;
                double terme_commun = 24.0 * (2 * inv_r12 - inv_r6) * inv_r2;

                double force_x = terme_commun * x_sous;
                double force_y = terme_commun * y_sous;

                tab_par[i].fx += force_x;
                tab_par[j].fx -= force_x;
                tab_par[i].fy += force_y;
                tab_par[j].fy -= force_y;

                viriel += terme_commun * rij2;
            }
        }
    }

    return viriel; 
}


//-----------------------------{linspace fonction}---------------------------------------

double* linspace(double x1, double x2, int n){
    double *x = (double *) calloc(n, sizeof(double)); 
    double step = (x2 - x1)/(double)(n-1);
    for(int i=0; i<n; i++){
        x[i]=x1+((double)i*step);
    }
    return x; 
}

//-----------------------------{arange fonction}---------------------------------------

// Fonction arange sans size en argument
double* arange(double start, double stop, double step) {
    // Calculer la taille du tableau
    int size = (int)((stop - start) / step);
    if (size <= 0) {
        return NULL;  // Retourne NULL si aucune valeur à générer
    }

    // Allouer de la mémoire pour le tableau
    double* array = (double*)malloc(size * sizeof(double));
    if (array == NULL) {
        printf("Erreur d'allocation mémoire.\n");
        exit(1);
    }

    // Remplir le tableau avec les valeurs
    for (int i = 0; i < size; i++) {
        array[i] = start + i * step;
    }

    return array;
}

//-----------------------------{main}----------------------------------------------------------------
int main() 
{
    // Initialiser le générateur de nombres aléatoires avec une graine pour avoir des résultats reproductibles
    srand(time(NULL));

    // Creation d'un tableau de particules
    Particule tab_particules[nbx_particules];
    Particule tab_particules_temporaire[nbx_particules]; 


    double pression;
    double viriel; 

    int taille = 999; 
    double* temperature_variations_rechauffement = linspace(0.01,10,taille); 
    double* temperature_variations_refroidissement = arange(10,0,-0.01); 


    // Initialisation des particules avec une configuration aleatoire
    initialiser_Configuration_Cristalline(tab_particules); 
    
    // Ouvre le fichier "positions_data.txt" en mode écriture ("w")
    // Ce fichier stockera les positions des particules au fil du temps dans le format suivant:
    // temps   x1 y1    x2 y2    x3 y3    x4 y4.....
    FILE *energie_file = fopen("energie_reversibilite_L40.txt", "w"); 

    int tentative = 0;
    int tentative_acceptee=1; 
    int index_modidi=0; 
    int cycle_externe=1; 
    int cycle_externe_refroidissement=1; 


    const int nombre_cycle_total = 120000+taille+20000+taille+120000;
    const int nombre_cycle_post_chauffe = 120000+taille; 
    const int nombre_cycle_post_chauffe_equilibre = 120000+taille+20000;
    const int nombre_cycle_post_refroidissement = 120000+taille+20000+taille; 

    double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
    enregistrer_energie(energie_file, energie_cycle, 0); 

    for (int cycle=0 ; cycle<nombre_cycle_total; cycle++){


            if (cycle <= 120000){

                printf("phase 1, équilibre préchauffe : %d/120000\n", cycle); 

                for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules_temporaire[i]=tab_particules[i];
                    }

                    int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                    double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi, 0.01); 
                    //printf("%f\n", proba_finale);

                    double x= ((double)rand() / RAND_MAX);

                    if (proba_finale >= x){

                        tentative_acceptee += 1; 

                        for (int i=0; i<nbx_particules; i++){
                            tab_particules[i]=tab_particules_temporaire[i];
                        }
                    }

                    tentative += 1; 
                }

                double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
                enregistrer_energie(energie_file, energie_cycle, cycle); 

            } else if (120000 < cycle && cycle <= nombre_cycle_post_chauffe){

                printf("phase 2, chauffe : %d/121000\n", cycle);  

                for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules_temporaire[i]=tab_particules[i];
                    }

                    int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                    double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi, temperature_variations_rechauffement[cycle_externe]); 

                    double x= ((double)rand() / RAND_MAX);

                    if (proba_finale >= x){

                        tentative_acceptee += 1; 

                        for (int i=0; i<nbx_particules; i++){
                            tab_particules[i]=tab_particules_temporaire[i];
                        }
                    }

                    tentative += 1; 
                }

                double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
                enregistrer_energie(energie_file, energie_cycle, cycle);

                cycle_externe +=1; 
            } else if (nombre_cycle_post_chauffe < cycle && cycle <= nombre_cycle_post_chauffe_equilibre){

                printf("phase 3, post chauffe equilibre : %d/141000\n", cycle);  

                for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules_temporaire[i]=tab_particules[i];
                    }

                    int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                    double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi, 10); 

                    double x= ((double)rand() / RAND_MAX);

                    if (proba_finale >= x){

                        tentative_acceptee += 1; 

                        for (int i=0; i<nbx_particules; i++){
                            tab_particules[i]=tab_particules_temporaire[i];
                        }
                    }

                    tentative += 1; 
                }

                double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
                enregistrer_energie(energie_file, energie_cycle, cycle);

            } else if (nombre_cycle_post_chauffe_equilibre < cycle && cycle <= nombre_cycle_post_refroidissement){

                printf("phase 4, refroidissement : %d/142000\n", cycle);  

                for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules_temporaire[i]=tab_particules[i];
                    }

                    int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                    double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi, temperature_variations_refroidissement[cycle_externe_refroidissement]); 

                    double x= ((double)rand() / RAND_MAX);

                    if (proba_finale >= x){

                        tentative_acceptee += 1; 

                        for (int i=0; i<nbx_particules; i++){
                            tab_particules[i]=tab_particules_temporaire[i];
                        }
                    }

                    tentative += 1; 
                }

                double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
                enregistrer_energie(energie_file, energie_cycle, cycle);

                cycle_externe_refroidissement +=1; 

            } else {

                printf("phase 5, refroidissement equilibre : %d/262000\n", cycle);  

                for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules_temporaire[i]=tab_particules[i];
                    }

                    int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                    double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi, 0.01); 

                    double x= ((double)rand() / RAND_MAX);

                    if (proba_finale >= x){

                        tentative_acceptee += 1; 

                        for (int i=0; i<nbx_particules; i++){
                            tab_particules[i]=tab_particules_temporaire[i];
                        }
                    }

                    tentative += 1; 
                }

                double energie_cycle = calcul_energie_potentielle_terme_commun(tab_particules,index_modidi)+U_tail+calcul_energie_potentielle_terme_modifie(tab_particules,index_modidi);
                enregistrer_energie(energie_file, energie_cycle, cycle);
            }

        }
    fclose(energie_file);  // De même pour le fichier positions_data.txt

    return 0;
 
}
