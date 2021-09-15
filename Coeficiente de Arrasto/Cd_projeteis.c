 /* Calculador de Coeficientes de Arrasto para projéteis subsônicos.
    Copyright (C) 2021  Mario Pereira (mv-pereira).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    */

/********************************************************
 * Programa para cálculo do Coeficiente de Arrasto      *
 * Cd (adimensional) para projéteis de dimensões        *
 * conhecidas utilizando parâmetros fornecidos pelos    *
 * fabricantes de projétis.                             *
 *                                                      *
 * É nessário saber a perda de velocidade do projétil   *
 * em determinada distância.                            *
 *                                                      *
 * Mario Pereira                                        *
 * Perito Criminal - DEPP/ICPAS - Pernambuco            *
 ********************************************************/

#include <stdio.h>
#include <math.h>
#include <locale.h> //Utilizando caracteres e acentuação da língua portuguesa.
#include <stdlib.h> //Para função exit() na condicional da abertura do arquivo;


#define OMEGA 0.000072921   //Taxa de rotação da terra em "rad/s".
#define PARADA 0.0174533    // Critério de parada para ajuste de angulação. 0.0174533 rad = 1º.
#define H 0.0001            //passo da iteração do Runge-Kutta.
#define DEBUG 0             //DEBUG == 1 permite iniciar o programa com parâmetros pré setados.

//Estrutura do Projétil

enum sentido_rotacao {Dextrogiro,Levogiro};

struct caracteristicas_do_projetil { //Todas medidas com base no projétil
    enum sentido_rotacao rotacao;
    double massa;
    double diametro;
    double coef_arrasto;
};

struct prjt {
    double x, y, z;
    double vx, vy, vz;
    double taxa_de_subida, rumo; //inclinacao e inclinação lateral instantânea
    double latitude, longitude;
    struct caracteristicas_do_projetil propriedades;
};

//Estrutura do Vento
struct vento {
    double velocidade, direcao;
    double x,y,z;           //Componentes x,y,z na direção de deslocamento principal (downrange) do projétil.
    double norte,leste;     //Componentes Norte e Leste da direção do vento.
};

//Estrutura da Edificação
struct edificacao {
    double latitude;
    double longitude;
    double altura;
};

//Estrutura do Disparo
enum origem_disparo {Nivel_do_Mar, Edificacao};

struct disparo {
    enum origem_disparo origem;
    double latitude;
    double longitude;
    double altura;
    double velocidade;      //Estimado pela característica do projétil.
    double azimute;
    double theta;
};

//Estrutura da Impactação - Todas medidas na cena de crime.
struct impactacao {
    double latitude;
    double longitude;
    double altura;
    double phi;
    double azimute;
};

/****************************************************************************
 *Aproximacao exponencial para Densidade                                    *
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


/****************************************************
 * Função Principal: Requer alguns iniciadores para *
 *                   calcular Cd, o Coeficiente     *
 *                   de Arrasto.                    *
 *                                                  *
 ****************************************************/

int main(){
    setlocale(LC_ALL, "Portuguese"); //Utilizando caracteres e acentuação da língua portuguesa.

    printf("\nPrograma para cálculo do Coeficiente de Arrasto Cd (adimensional) para projéteis de dimensões conhecidas utilizando parâmetros fornecidos pelos fabricantes de projétis.\n");
    
    double t, g;                       //Tempo e Aceleração da gravidade. t1 representa o tempo no momento n+1.
    double kappa, kappaSdensidade;     //Valores da característica relacionada ao arrasto.

    double distancia,velocidadeF,delta_v;
    
    struct prjt projetil;               //Estrutura do projétil.
    struct vento w;                     //Definição da struct do vento.
    struct impactacao impacto;          //Estrutura da Impactação.
    struct disparo tiro;                //Estrutura do Tiro.

    tiro.latitude = -23.698308;         //Latitude da Fábrica da CBC em São Paulo
    
#if DEBUG   //Valores padrão de .40 ETPP
tiro.velocidade = 302.0;
distancia = 100.0;
velocidadeF = 274.0;
projetil.propriedades.massa = 11.66/1000;
projetil.propriedades.diametro = M_PI*powl((10/1000.0),2)/4;
printf ("\n*\t*\tDEBUG Ativado.\t*\t*\nValores prefixados para um projétil calibre .40 ETPP.\n");

#else
/********************************
 * Características do projetil  *
 ********************************/

    printf("\nDigite a velocidade inicial (em m/s) do projétil:\n");
    scanf("%lf", &tiro.velocidade);

    printf("\nDigite para qual distancia (em m) o Cd será calculado:\n");
    scanf("%lf", &distancia);
    
    printf("\nDigite a velocidade final (em m/s) do projétil:\n");
    scanf("%lf", &velocidadeF);

    printf("\nDigite a massa (em g) do projétil:\n");
    scanf("%lf", &projetil.propriedades.massa);
    projetil.propriedades.massa = projetil.propriedades.massa/1000.0; //No SI, massa em Kg

    printf("\nDigite o diâmetro (em mm) do projétil Disparado:\n");
    scanf("%lf", &projetil.propriedades.diametro);
    // A variável chama-se a para cálculo da área transversal.
    projetil.propriedades.diametro = M_PI*powl((projetil.propriedades.diametro/1000.0),2)/4;

#endif
    
    //Valor médio estimado para início dos cálculos.
    projetil.propriedades.coef_arrasto = 0.2;

    //Aceleração da gravidade na latitude. (em m/s^2)
    g = 9.780327*(1+0.0053024*sin(tiro.latitude)*sin(tiro.latitude) - 0.0000058*sin(2*tiro.latitude)*sin(2*tiro.latitude));

/********************************
 * Ponto de Partida do GOTO após*
 * correção do c = Constante de *
 * Arrasto                      *
 ********************************/

    A:
    
    kappaSdensidade = (projetil.propriedades.coef_arrasto/projetil.propriedades.massa)*(projetil.propriedades.diametro/2); //Ordem alterada de c*a/(2.0*m) para evitar Underflow.
    //kappaSdensidade -> κ/densidade_ar. Note que a massa já entra na constante. Densidade do ar será calculado dentro do laço...
    
/********************************************************
 * Condições Iniciais: Dados dispostos considerando um  *
 *                     disparo a um metro de altura     *
 *                     a 0 grau (projetil.vy=0).        *
 *                                                      *
 ********************************************************/
    t=0.0;
    projetil.x=0.0;
    projetil.y=1;
    projetil.z=0.0;
    projetil.vx=tiro.velocidade*cos(0);
    projetil.vy=tiro.velocidade*sin(0);
    projetil.vz=0.0;                    //Downrange continua sendo no eixo x;

    projetil.taxa_de_subida = atan2 (projetil.vy,projetil.vx);
    projetil.rumo = tiro.azimute + atan2 (projetil.vz,projetil.vx);

    w.x=0.0; //Condições iniciais do vento.
    w.y=0.0;
    w.z=0.0;
    
/************************************
 * Início do laço para cálculo RK4  *
 *                                  *
 ************************************/
    
    while (projetil.x < distancia){
        t += H;
        projetil.x += runge_kutta(&pos_x, &projetil, &w, H, 0, 0)*H;
        projetil.y += runge_kutta(&pos_y, &projetil, &w, H, 0, 0)*H;
        projetil.z += runge_kutta(&pos_z, &projetil, &w, H, 0, 0)*H;
        
        projetil.rumo = tiro.azimute + atan2 (projetil.vz,projetil.vx);

        kappa = kappaSdensidade*densidade_ar(projetil.y);

        // Nos cálculos iterativos da velocidade nos eixos, as variáveis: w.x,w.y e w.z (Velocidade do vento), latitude e azimute foram desconsideradas.

        projetil.vx += runge_kutta(&w_vx, &projetil, &w, H, kappa, g)*H;
        projetil.vy += runge_kutta(&w_vy, &projetil, &w, H, kappa, g)*H;
        projetil.vz += runge_kutta(&w_vz, &projetil, &w, H, kappa, g)*H;

        projetil.taxa_de_subida = atan2 (projetil.vy,projetil.vx);

        }

/************************************
 * Etapa de correção do Coeficiente *
 * de atrito inicial.               *
 *                                  *
 ************************************/    

    delta_v = projetil.vx - velocidadeF; // Δv precisa ser menor do que uma quantidade ε. Aqui escolhido: 0.01 m/s.

    if ( fabs (delta_v) > 0.01){

    if ( delta_v > 0 ){               //Projétil terminou com mais velocidade que a Vf dada pela fabricante, após x metros de Downrange.
        projetil.propriedades.coef_arrasto += 0.00001;               //Adicione coeficiente de Arrasto.
        goto A;
        }
    else{                             //Projétil terminou com menos velocidade que a Vf dada pela fabricante, após x metros de Downrange.
        projetil.propriedades.coef_arrasto -= 0.00001;               //Retire coeficiente de Arrasto.
        goto A;
        }
    }
    
    printf ("\nConsiderando uma perda de velocidade de %.2f m/s em %.2f m,\nO coeficiente de Arrasto vale: Cd=%f.\n",tiro.velocidade-velocidadeF,distancia,projetil.propriedades.coef_arrasto);
    
    return 0;
}
