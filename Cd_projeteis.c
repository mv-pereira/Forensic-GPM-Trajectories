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

struct prjt { //estrutura que guarda a posição tridimensional e sua velocidade
    double x,y,z,vx,vy,vz,latitude,azimute;
};

//Estrutura do Vento
struct vento { //Estrutura para guardar a velocidade do vento nos eixos.
    double x,y,z;
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

/************************************************
 * Funções para cálculo de Runge-Kutta 4a Ordem *
 ************************************************/

double pos(double v_na_direcao, double h){
    double k,kx1,kx2,kx3,kx4;
    kx1 = v_na_direcao;
    kx2 = v_na_direcao+kx1*(h/2);
    kx3 = v_na_direcao+kx2*(h/2);
    kx4 = v_na_direcao+kx3*h;
    k = (1/6.0)*(kx1 + 2*kx2 + 2*kx3 + kx4);
    return k;
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

double w_vx (double k, struct prjt *projetil, double correcao, struct vento *w){
    return  (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vx +correcao -w->x) + 2*OMEGA*( -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(sin (projetil->azimute)) -(projetil->vz +correcao -w->z)*(sin (projetil->latitude))));
}

double kvx (struct prjt *projetil, struct vento *w, double h, double kappa){
    double k,kvx1,kvx2,kvx3,kvx4;
    kvx1 = w_vx(kappa, projetil, 0, w);
    kvx2 = w_vx(kappa, projetil, kvx1*(h/2), w);
    kvx3 = w_vx(kappa, projetil, kvx2*(h/2), w);
    kvx4 = w_vx(kappa, projetil, kvx3*(h/2), w);
    k = (1/6.0)*(kvx1 + 2*kvx2 + 2*kvx3 + kvx4);
    return k;
}

/*****************************************************************
 * funcao Auxiliar w_vy para Vel y (Altura): Esta função apenas  *
 *                                           calculra o 'k' para *
 *                                           projetil_1.vy.      *
 *                                                               *
 *****************************************************************/


double w_vyy (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vy +correcao -w->y) -g + 2*OMEGA*( (projetil->vx +correcao -w->x)*(cos (projetil->latitude))*(sin (projetil->azimute)) +(projetil->vz +correcao -w->z)*(cos (projetil->latitude))*(sin (projetil->azimute)) ) );
} 


double kvy (struct prjt *projetil, struct vento *w, double h, double kappa, double g){
    double k,kvy1,kvy2,kvy3,kvy4;
    kvy1 = w_vyy(kappa, projetil, 0, w,g);
    kvy2 = w_vyy(kappa, projetil, kvy1*(h/2), w,g);
    kvy3 = w_vyy(kappa, projetil, kvy2*(h/2), w,g);
    kvy4 = w_vyy(kappa, projetil, kvy3*(h/2), w,g);
    k = (1/6.0)*(kvy1 + 2*kvy2 + 2*kvy3 + kvy4);
    return k;
}


/******************************************************************
 * funcao Auxiliar para Vel z (Deriva/drift): Esta função apenas  *
 *                                            calculra o 'k' para *
 *                                            vz1.                *
 *                                                                *
 ******************************************************************/

double w_vz (double k, struct prjt *projetil, double correcao, struct vento *w){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vz +correcao -w->z) + 2*OMEGA*( (projetil->vx +correcao -w->x)*(sin (projetil->latitude)) -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(cos (projetil->azimute)) )   );
}


double kvz (struct prjt *projetil, struct vento *w, double h, double kappa){
    double k,kvz1,kvz2,kvz3,kvz4;
    kvz1 = w_vz(kappa, projetil, 0, w);
    kvz2 = w_vz(kappa, projetil, kvz1*(h/2), w);
    kvz3 = w_vz(kappa, projetil, kvz2*(h/2), w);
    kvz4 = w_vz(kappa, projetil, kvz3*(h/2), w);
    k = (1/6.0)*(kvz1 + 2*kvz2 + 2*kvz3 + kvz4);
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
    
    double v0;                                      //Velocidade Inicial
    double t,t1,g;                                  //Tempo e Aceleração da gravidade. t1 representa o tempo no momento n+1.
    double c, a, massa, kappa, kappaSdensidade;     //Valores de entrada (Coef. Arrasto; densidade do Ar; área de sec trasnversal do objeto, massa, constante).
    double latitude = -23.698308;                   //Latitude da Fábrica da CBC em São Paulo
    double distancia,velocidadeF,delta_v;
    
    struct prjt projetil, projetil_1;   //projetil_1 é o projétil para n+1 | Definição dessa struct no começo do programa.
    struct vento w;                     //Definição da struct do vento.

    
#if DEBUG   //Valores padrão de .40 ETPP
v0 = 302.0;
distancia = 100.0;
velocidadeF = 274.0;
massa = 11.66/1000;
a = M_PI*powl((10/1000.0),2)/4;
printf ("\n*\t*\tDEBUG Ativado.\t*\t*\nValores prefixados para um projétil calibre .40 ETPP.\n");

#else
/********************************
 * Características do projetil  *
 ********************************/

    printf("\nDigite a velocidade inicial (em m/s) do projétil:\n");
    scanf("%lf", &v0);

    printf("\nDigite para qual distancia (em m) o Cd será calculado:\n");
    scanf("%lf", &distancia);
    
    printf("\nDigite a velocidade final (em m/s) do projétil:\n");
    scanf("%lf", &velocidadeF);

    printf("\nDigite a massa (em g) do projétil:\n");
    scanf("%lf", &massa);
    massa = massa/1000.0; //No SI, massa em Kg

    printf("\nDigite o diâmetro (em mm) do projétil Disparado:\n");
    scanf("%lf", &a);
    // A variável chama-se a para cálculo da área transversal.
    a = M_PI*powl((a/1000.0),2)/4;

#endif
    
    //Valor médio estimado para início dos cálculos.
    c = 0.2;

    //Aceleração da gravidade na latitude. (em m/s^2)
    g = 9.780327*(1+0.0053024*sin(latitude)*sin(latitude) - 0.0000058*sin(2*latitude)*sin(2*latitude));

/********************************
 * Ponto de Partida do GOTO após*
 * correção do c = Constante de *
 * Arrasto                      *
 ********************************/

    A:
    
    kappaSdensidade = c*a/(2.0*massa);
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
    projetil.vx=v0*cos(0);
    projetil.vy=v0*sin(0);
    projetil.vz=0.0;                    //Downrange continua sendo no eixo x;
    projetil.latitude = -23.698308;     //Latitude da Fábrica da CBC em São Paulo
    projetil.azimute = 0;
    projetil_1 = projetil;
    
    w.x=0.0; //Condições iniciais do vento.
    w.y=0.0;
    w.z=0.0;
    
/************************************
 * Início do laço para cálculo RK4  *
 *                                  *
 ************************************/
    
    while (projetil_1.x < distancia){
        t1 = t + H;
        projetil_1.x = projetil.x + pos(projetil.vx,H)*H;
        projetil_1.y = projetil.y + pos(projetil.vy,H)*H;
        projetil_1.z = projetil.z + pos(projetil.vz,H)*H;
        
        kappa = kappaSdensidade*densidade_ar(projetil.y);

        // Nos cálculos iterativos da velocidade nos eixos, as variáveis: w.x,w.y e w.z (Velocidade do vento), latitude e azimute foram desconsideradas.

        projetil_1.vx = projetil.vx + kvx(&projetil, &w, H, kappa)*H;        
        projetil_1.vy = projetil.vy + kvy(&projetil, &w, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + kvz(&projetil, &w, H, kappa)*H;
        
        // Atualização das variáveis.
        t = t1;
        projetil = projetil_1;
        }

/************************************
 * Etapa de correção do Coeficiente *
 * de atrito inicial.               *
 *                                  *
 ************************************/    

    delta_v = projetil_1.vx - velocidadeF; // Δv precisa ser menor do que uma quantidade ε. Aqui escolhido: 0.01 m/s.

    if ( fabs (delta_v) > 0.01){

    if ( delta_v > 0 ){               //Projétil terminou com mais velocidade que a Vf dada pela fabricante, após x metros de Downrange.
        c = c + 0.0001;               //Adicione coeficiente de Arrasto.
        goto A;
        }
    else{                             //Projétil terminou com menos velocidade que a Vf dada pela fabricante, após x metros de Downrange.
        c = c - 0.0001;               //Retire coeficiente de Arrasto.
        goto A;
        }
    }
    
    printf ("\nConsiderando uma perda de velocidade de %.2f m/s em %.2f m,\nO coeficiente de Arrasto vale: Cd=%f.\n",v0-velocidadeF,distancia,c);
    
    return 0;
}
