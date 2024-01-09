#include "projetil.h"



#define OMEGA 0.000072921   //Taxa de rotação da terra em "rad/s".



/****************************************************************************
 *Aproximacao exponencial para Densidade do ar.                             *
 *https://en.wikipedia.org/wiki/Density_of_air#Exponential_approximation    *
 ****************************************************************************/

double densidade_ar (double altura){
    return 1.22501236323*exp(-altura/10400);
}

/****************************************
 * Funções trigonométricas auxiliares   *
 ****************************************/

double sec(double alpha){
    return 1/(cos (alpha));
}
double arcsec(double x){
    return acos(1/x); //arcsec t = arccos(1/t).
}
/************************************************************
 * Fórmula de Haversine                                     *
 * Importante equação usada em navegação, fornecendo        *
 * distâncias entre dois pontos de uma esfera a partir      *
 * de suas latitudes e longitudes (Recebe em grau, mas      *
 * precisa converter para RADIANOS).                        *
 * https://pt.wikipedia.org/wiki/F%C3%B3rmula_de_haversine  *
 ************************************************************/
double raioTerraLatitude (double latitude){ //https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius
    double a, b;
    a = 6378137.0; //Raio Equatorial em m.
    b = 6356752.3; //Raio Polar em m.
    return sqrt( (pow(a*a*cos(latitude),2) + pow(b*b*sin(latitude),2) ) / (pow(a*cos(latitude),2) + pow(b*sin(latitude),2) ) );
}

double haversine(double lat1, double long1, double lat2, double long2){
    double raio;
    raio = raioTerraLatitude ((lat1+lat2)/2); //Raio da terra.
    lat1= M_PI*lat1/180;
    lat2= M_PI*lat2/180;
    long1= M_PI*long1/180;
    long2= M_PI*long2/180;
    return 2*raio*asin(sqrt( pow( sin((lat2-lat1)/2),2) + cos(lat1)*cos(lat2)*pow(sin((long2-long1)/2),2)));
}


/********************************************************************
 * Cáluclo da distância percorrida pelo projétil convertida         *
 * para graus latitude/longitude tendo em vista sua localização     *
 * relativa às edificações e impactações.                           *
 *                                                                  *
 * "distLatGraus" é a "distância" em graus entre a latitude do      *
 * disparo e a latitude do ponto solicitado.                        *
 * "distLongGraus" é a "distância" em graus entre a longitude do    *
 * disparo e a longitude do ponto solicitado.                       *
 *                                                                  *
 * Constante "0.000008983152098" tem unidades de º/m.               *
 ********************************************************************/

double distLatGraus (const struct prjt *projetil, const struct disparo *tiro){
    return 0.000008983152098*(projetil->x*cos(tiro->azimute) - projetil->z*sin(tiro->azimute));
}

double distLongGraus (const struct prjt *projetil, const struct impactacao *impacto, const struct disparo *tiro){
    return 0.000008983152098*cos(impacto->latitude)*(projetil->x*sin(tiro->azimute) + projetil->z*cos(tiro->azimute) );
}


/****************************************************************************
 *                        Spindrift aproximado                              *
 * Valores estimados para projéteis subsônicos a partir de                  *
 * https://theoverwatch.wixsite.com/theoverwatch/post/spin-drift            *
 * sg = 2.42563 para projétil subsônico (300 Blackout using subsonic ammo)  *
 * sg = 1.52551 para projétil supersônico .308                              *
 ****************************************************************************/
double spindrift (double tempo){
    double sg,drift;
    sg = 2.42563;
    drift = 1.25*(sg+1.2)*pow(tempo,1.83);  //Resultado em polegada.
    return 2.54*drift;                      //Resultado em centímetro.
}


/****************************************************************
 * funções auxiliares para cálculo de posição no RK.            *
 ****************************************************************/

double pos_x (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return  projetil->vx;
}

double pos_y (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return  projetil->vy;
}

double pos_z (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return  projetil->vz;
}

/********************************************************************
 * funcao Auxiliar w_vx para Vel x (Downrange): Esta função apenas  *
 *                                              calculra o 'k' para *
 *                                              projetil_1.vx.      *
 * Esta função w_vx (y ou z) recebe o endereço da estrutura projétil*
 * a correção da iteração kvx1*(h/2) e o endereço da estrutura      *
 * da velocidade do vento                                           *
 *                                                                  *
 ********************************************************************/
double w_vx (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vx +correcao -w->x) + 2*OMEGA*( -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(sin (projetil->rumo)) -(projetil->vz +correcao -w->z)*(sin (projetil->latitude))));
}

/*****************************************************************
 * funcao Auxiliar w_vy para Vel y (Altura): Esta função apenas  *
 *                                           calculra o 'k' para *
 *                                           projetil_1.vy.      *
 *                                                               *
 *****************************************************************/
double w_vy (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vy +correcao -w->y) -g + 2*OMEGA*( (projetil->vx +correcao -w->x)*(cos (projetil->latitude))*(sin (projetil->rumo)) +(projetil->vz +correcao -w->z)*(cos (projetil->latitude))*(sin (projetil->rumo)) ) );
}

/******************************************************************
 * funcao Auxiliar para Vel z (Deriva/drift): Esta função apenas  *
 *                                            calculra o 'k' para *
 *                                            projetil_1.vz1.     *
 *                                                                *
 ******************************************************************/
double w_vz (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vz +correcao -w->z) + 2*OMEGA*( (projetil->vx +correcao -w->x)*(sin (projetil->latitude)) -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(cos (projetil->rumo)) )   );
}

/************************************************
 * Função para cálculo de Runge-Kutta 4a Ordem  *
 ************************************************/
double runge_kutta (double (*funcao) (double, struct prjt (*), double, struct vento (*), double), struct prjt *projetil, struct vento *w, double passo, double kappa, double g){
    double k,k1,k2,k3,k4;
    
    k1 = funcao(kappa, projetil, 0, w, g);
    k2 = funcao(kappa, projetil, k1*(passo/2), w, g);
    k3 = funcao(kappa, projetil, k2*(passo/2), w, g);
    k4 = funcao(kappa, projetil, k3*(passo/2), w, g);
    k = (1/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    return k;
}


/********************************************************
 * Ajuste de θ inicial para coincidência de ϕ (final)   *
 *                                                      *
 ********************************************************/

 double ajustar_theta (double phi_final, double phi_medido, double theta){
    double grau_sobre_cem_rad = 0.0001745329;
    static int escala = 2;
    double diferenca_escalada = fabs (phi_final - phi_medido)/escala;

    if ( fabs (phi_medido) < 0.3) grau_sobre_cem_rad = 0.0000174532; //a correção passa a ser de 0.001 grau
    if (phi_medido == 0) grau_sobre_cem_rad = 0.0000017453;          //Necessário para aumento de precisão se phi == 0;
    if ( phi_medido < phi_final ){              //Projétil terminou com angulação maior que o φ medido, indicando que o disparo foi mais baixo.
        if (phi_medido >= 0) theta -= diferenca_escalada + grau_sobre_cem_rad;
        else                 theta += diferenca_escalada + grau_sobre_cem_rad;
    }
    else{                                   //Projétil terminou com angulação menos que o φ medido, indicando que o disparo foi mais alto.
        /*phi_final < phi_medido*/
        if (phi_medido >= 0) theta += diferenca_escalada + grau_sobre_cem_rad;
        else                 theta -= diferenca_escalada + grau_sobre_cem_rad;
        }
    escala++; // diferença escalada a cada iteração diminui a contribuição (φf-φm)/2, depois /3, /4...

    return theta;
 }

 /*******************************************************
 * Ajuste de AZ inicial para coincidência de impacto.azimute (γ)  *
 * γ = azimute_medido                                   *
 ********************************************************/

double ajuste_AZ(double azimute_disparo, double azimute_final, double azimute_medido){

    double grau_sobre_cem_rad = 0.0001745329;
    static int escala = 2;
    double diferenca_escalada = fabs (azimute_final - azimute_disparo)/escala;

    if (azimute_final > azimute_medido) {
        azimute_disparo -= diferenca_escalada + grau_sobre_cem_rad;
    } else {
        azimute_disparo += diferenca_escalada + grau_sobre_cem_rad;
    }
    escala++; // diferença escalada a cada iteração diminui a contribuição (azf-azm)/2, depois /3, /4...

    return azimute_disparo;
}