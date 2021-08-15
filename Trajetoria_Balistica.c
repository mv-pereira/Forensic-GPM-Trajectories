 /* Calculador de Trajetórias para projéteis subsônicos.
  
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
 * Programa para cálculo de Trajetória Balística        *
 * utilizando método de massa pontual, considerando:    *
 * - Variação de densidade com altura;                  *
 * - Variação de aceleração de gravidade com latitude;  *
 * - Efeito Coriolis;                                   *
 * - Efeito da deriva ocasionado pelo vento;            *
 *                                                      *
 * Mario Pereira                                        *
 * Perito Criminal - DEPP/ICPAS - Pernambuco            *
 ********************************************************/


#include <stdio.h>
#include <math.h>
#include <locale.h> //Utilizando caracteres e acentuação da língua portuguesa.
#include <stdlib.h> //Para função exit() na condicional da abertura do arquivo;
#include <stdbool.h>
#include <time.h>   //Para velocidade de funções - DEBUG


#define OMEGA 0.000072921   //Taxa de rotação da terra em "rad/s".
#define PARADA 0.00174533    // Critério de parada para ajuste de angulação. 0.00174533 rad = 0.1º.
#define H 0.0001            //passo da iteração do Runge-Kutta.
#define DEBUG 1

//Estrutura do Projétil

struct prjt { //estrutura que guarda a posição tridimensional e suas velocidades
    double x,y,z,vx,vy,vz,latitude,azimute;
};

//Estrutura do Vento
struct vento { //Estrutura para guardar a velocidade do vento nos eixos.
    double x,y,z;
};

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

double distLatGraus (struct prjt *projetil){
    return 0.000008983152098*(projetil->x*cos(projetil->azimute) - projetil->z*sin(projetil->azimute));
}

double distLongGraus (struct prjt *projetil, double latitude){
    return 0.000008983152098*cos(latitude)*(projetil->x*sin(projetil->azimute) + projetil->z*cos(projetil->azimute) );
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
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vx +correcao -w->x) + 2*OMEGA*( -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(sin (projetil->azimute)) -(projetil->vz +correcao -w->z)*(sin (projetil->latitude))));
}

/*****************************************************************
 * funcao Auxiliar w_vy para Vel y (Altura): Esta função apenas  *
 *                                           calculra o 'k' para *
 *                                           projetil_1.vy.      *
 *                                                               *
 *****************************************************************/
double w_vy (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vy +correcao -w->y) -g + 2*OMEGA*( (projetil->vx +correcao -w->x)*(cos (projetil->latitude))*(sin (projetil->azimute)) +(projetil->vz +correcao -w->z)*(cos (projetil->latitude))*(sin (projetil->azimute)) ) );
}

/******************************************************************
 * funcao Auxiliar para Vel z (Deriva/drift): Esta função apenas  *
 *                                            calculra o 'k' para *
 *                                            projetil_1.vz1.     *
 *                                                                *
 ******************************************************************/
double w_vz (double k, struct prjt *projetil, double correcao, struct vento *w, double g){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vz +correcao -w->z) + 2*OMEGA*( (projetil->vx +correcao -w->x)*(sin (projetil->latitude)) -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(cos (projetil->azimute)) )   );
}

/************************************************
 * Função para cálculo de Runge-Kutta 4a Ordem  *
 ************************************************/
double runge_kutta (double (*funcao) (double, struct prjt (*), double, struct vento (*), double), struct prjt *projetil, struct vento *w, double passo, double kappa, double g){
    double k,k1,k2,k3,k4;
    double inclinacao_lateral;
    inclinacao_lateral = atan (projetil->vz/projetil->vx);
    projetil->azimute += inclinacao_lateral;

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

    if ( fabs (phi_medido) < 0.3) grau_sobre_cem_rad = grau_sobre_cem_rad/10; //a correção passa a ser de 0.001 grau

    if ( phi_medido < phi_final ){              //Projétil terminou com angulação maior que o φ medido, indicando que o disparo foi mais baixo.
        if (phi_medido >= 0) theta -= grau_sobre_cem_rad;
        else                 theta += grau_sobre_cem_rad;
    }
    else{                                   //Projétil terminou com angulação menos que o φ medido, indicando que o disparo foi mais alto.
        /*phi_final < phi_medido*/
        if (phi_medido >= 0) theta += grau_sobre_cem_rad;
        else                 theta -= grau_sobre_cem_rad;
        }

    return theta;
 }

 /*******************************************************
 * Ajuste de AZ inicial para coincidência de gamma (γ)  *
 * γ = azimute_medido                                   *
 ********************************************************/

double ajuste_AZ(double azimute_disparo, double inclinacao_lateral, double azimute_medido){

    double grau_sobre_cem_rad = 0.0001745329;

    double azimute_final = azimute_disparo + inclinacao_lateral;

    if (azimute_final > azimute_medido) {
        azimute_disparo -= grau_sobre_cem_rad;
    } else {
        azimute_disparo += grau_sobre_cem_rad;
    }

    return azimute_disparo;
}



/****************************************************
 * Função Principal: Trajetória Balística           *
 *                                                  *
 ****************************************************/

int main(){
    setlocale(LC_ALL, "Portuguese"); //Utilizando caracteres e acentuação da língua portuguesa.

    FILE *arquivo;
    arquivo = fopen("data","w");
    
#if DEBUG
    /*DEBUG para geração de arquivos com informações pertinentes.*/
        FILE *debug;
        debug = fopen("debug","w");
    /*Cálculo de duração de funções*/
    clock_t ciclos_cpu;
#endif

    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        exit(1);
    }
    
    int n;                              //Contador.
    int dextrogiro;                     //variável "booleana" para critério do laço, e se o projétil é dextrogiro ou levogiro.
    double altura, phi, v;              //Valores de entrada.
    double xs[2];                       //x+ e x- Raízes do sistema sem arrasto ("s" sem arrasto).
    double theta;                       //Angulo de disparo.
    double t,t1,t_total,g;              //Tempo, Tempo em n+1 e aceleração da gravidade.
    double wx, wy, wz,vento,ang;        //Variáveis relacionadas à velocidade do vento e sua direção.
    double vN, vE;                      //Variáveis que participarão da decomposição do vento nos eixos x e z.
    double latitude, longitude, azimute;            //Variáveis realcionadas ao sistema de massa pontual.
    double c, ro, a, m, kappa, kappaSdensidade;     //Valores de entrada (Coef. Arrasto; densidade do Ar; área de sec trasnversal do objeto, massa, constante).
    double inclinacao, delta_phi;                                   //parâmetro de inclinação do projétil com arrasto e sua diferença com a medida.
    double inclinacao_lateral, delta_inclinacao_lateral, gamma;     //parâmetro para ajuste da inclinação lateral.
    double origemNMM[2];    //Vetor Origem (latitude e Longitude ao Nível Médio do Mar (NMM) do disparo (se ele fosse disparado ao NMM). Será gravado após o primeiro conjunto de iterações.
    
/*****************************************************************
 * Declaração das variaveis da possivel Edificação no meio do    *
 * caminho. Inicia-se o valor sem edificação e depois corrige    *
 * haviaEdf == 0 -> Não havia                                    *
 * haviaEdf == 1 -> Havia.                                       *
 *****************************************************************/
    double latitudeEdf,longitudeEdf,alturaEdf;
    int haviaEdf = 0;
    double altura_disparo = 0; //Altura do disparo inicia como sendo 0, ao NMM, depois corrige caso haja edificações na trajetória do projétil.
    double delta_y; //parâmetro para comparação entre a altura após atingir e a altura calculada após as iterações ao sair da edificação.
    double distanciaPredio_Impacataco; //Distancia entre a edificacao e a impactacao a ser calculada caso haja edificações.
    double velocidade_final;           //Velocidade do projétil na impactactação.
    double latitude_disparo, longitude_disparo; //Em graus
    
    struct prjt projetil, projetil_1; //projetil_1 é o projétil para n+1 | Definição dessa struct no começo do programa.
    struct vento w;                     //Definição da struct do vento.

#if DEBUG //Valores padrão apenas para comparação

altura = 89;
phi = 4*M_PI/180;
gamma = 183*M_PI/180;
v = 230;
m=10.240/1000.0;
a=M_PI*powl((8.82/1000.0),2)/4; // A divisão por quatro leva em conta o raio.
dextrogiro=1;
c = 0.235800;
vento = 10/3.6;
ang = (100+180.0)*M_PI/180.0;
latitude = -8.127727*M_PI/180.0;
g = 9.780327*(1+0.0053024*sin(latitude)*sin(latitude) - 0.0000058*sin(2*latitude)*sin(2*latitude)); //Açeleração da gravidade na latitude. (em m/s^2)
longitude = -34.898383;
azimute = gamma;

printf ("\n*\t*\tDEBUG Ativado.\t*\t*\n\t\tValores prefixados.\n\nPara sair da função DEBUG, mudar a definição de DEBUG para 0 no cabeçalho do programa e recompilar.\n\n");
//testando o prédio perto


#else

/********************************
 * Características do projetil  *
 ********************************/
    
    printf("Digite a altura (em metros) da impactação em relação ao disparo: ");
    scanf("%lf", &altura);

    printf("\nDigite o ângulo ϕ com a horizontal (em °):");
    printf("\nϕ ≈ 0 causa um bug conhecido. Se a impactação for descendente, ϕ < 0: ");
    scanf("%lf", &phi);
    phi = phi*M_PI/180; //Conversão de grau para radianos.
    
    printf("\nDigite o ângulo γ com o Norte (em °) (azimute): ");
    scanf("%lf", &gamma);
    gamma = gamma*M_PI/180; //Conversão de grau para radianos.
    
    printf("\nDigite a velocidade inicial (em m/s) do projétil: ");
    scanf("%lf", &v);

    printf("\nDigite a massa do projétil (em g): ");
    scanf("%lf", &m);
    m=(1.0*m)/1000.0; //[m] = Kg no SI
    
    printf("\nDigite o diâmetro do projétil já disparado (em mm): ");
    scanf("%lf", &a);
    a=M_PI*powl((a/1000.0),2)/4; //"a" assume agora o valor da Área de seção transversal para cálculo de kappa.
    
    printf("\nDigite o Coeficiente de Arrasto Cd (aproximadamente 0.2 em casos subsônicos) do projétil: ");
    scanf("%lf", &c);
    
    printf("\nO projetil é dextrogiro ou levogiro? 1 - Dextrogiro.\t2 - Levogiro.\n");
    scanf("%d", &dextrogiro);
    (dextrogiro != 1) ? dextrogiro = 0 : 1;
 
    printf("\nDigite o módulo da velocidade do vento (em km/h): ");
    scanf("%lf", &vento);
    vento = vento/3.6; // Velocidade do vento precisa ser em m/s.

    printf("\nDigite a direção do vento em relação ao Norte (em °): ");
    scanf("%lf", &ang);
    ang = (ang+180.0)*M_PI/180.0; // A direção dos ventos é dada de onde ele vem, não para onde sopra.
    
    printf("\nDigite a latitude decimal da impactacao (em °), Valor precisa ser negativo, se o disparo ocorrer Hemisfério Sul: ");
    scanf("%lf", &latitude);
    latitude = latitude*M_PI/180.0; //grau para radianos
    g = 9.780327*(1+0.0053024*pow(sin(latitude),2) - 0.0000058*pow (sin(2*latitude),2)); //Aceleração da gravidade na latitude. (em m/s^2) nvl do mar. g calculado antes de a la

    printf("\nDigite a longitude decimal da impactacao (em °): ");
    scanf("%lf", &longitude);
    
    /* A melhor estimativa para o Azimute inicial é o próprio gamma. Em condições normais de vento, não tem como divergir muito do γ */
    azimute = gamma;


#endif    

/****************************************************************************
 * Velocidade do vento é dada a partir de onde ele sopra e será decomposto  *
 * primeiramente nos eixos Norte e Leste e, posteriormente, nos eixos x     *
 * (Downrange) e z (desvio lateral)                                         *
 ****************************************************************************/

    vN = vento*cos(ang);
    vE = vento*sin(ang);

/************************************
 * Início do cálculo SEM arrasto    *
 * para estimativa do θ inicial     *
 * a ser utilizado nas iterações    *
 ************************************/    
    
    theta = arcsec (sqrtl(v*v*sec(phi)*sec(phi)/(v*v-2*g*altura)));
    printf("\nO ângulo θ do início do disparo com a horizontal considerado a partir do solo e em um sistema sem arrasto ou outras correcoes vale: %.2lf°\n",theta*180/M_PI);
    //Exibição da variável no printf convertida para grau.

    xs[0] = (powl(cos(phi),2)*(v*v - 2*g*altura)*( -fabs(tan(phi)) + sqrtl(v*v*powl(sec(phi),2)/(v*v-2*g*altura) - 1) ) ) /g ;
    xs[1] = (powl(cos(phi),2)*(v*v - 2*g*altura)*( +fabs(tan(phi)) + sqrtl(v*v*powl(sec(phi),2)/(v*v-2*g*altura) - 1) ) ) /g ;

    printf("\nConsiderando um sistema sem arrasto: x- = %.3lf m e x+ = %.3lf m\n",xs[0],xs[1]);
    
/********************************
 * Fim do cálculo sem Arrasto   *
 * onde foi obtido um θ inicial *
 ********************************/

    //Constante de arrasto/ro. Note que a massa já entra na constante. Densidade do ar será calculado nas funcoes de vx,vy,vz...
    kappaSdensidade = (c/m)*(a/2); //Ordem alterada de c*a/(2.0*m) para evitar Underflow.

    projetil.latitude=latitude;
    projetil.azimute =azimute;
    n = 0;                          //Contador para registrar quantas vezes toda a trajetória foi calculada.

#if DEBUG   /*Cálculo de duração de funções*/
ciclos_cpu = clock();
#endif

/********************************************************
 * Condições Iniciais: Dados dispostos considerando um  *
 *                     disparo a 0 metros de altura     *
 *                     e vy=v*sin(theta).               *
 *                                                      *
 ********************************************************/    
    
    A: //Ponto de partida do goto após correção do ângulo θ = theta.
    B: //Ponto de partida para correção do azimute do disparo.
    
    //Condições iniciais para variáveis com arrasto:
    t=0.0;
    projetil.x=0.0;
    projetil.y=0.0;
    projetil.z=0.0;
    projetil.vx=v*cos(theta);
    projetil.vy=v*sin(theta);
    projetil.vz=0.0;   //Downrange continua sendo no eixo x;
    //projetil.azimute será atualizado nos laços.
    projetil_1 = projetil;
   
    //Recálculo do vento pela variação do azimute inicial.
    //Mudança de base de Norte-Leste para X e Z
    w.x = vN*cos(projetil.azimute) + vE*sin(projetil.azimute);
    w.y = 0.0;
    w.z = -vN*sin(projetil.azimute) + vE*cos(projetil.azimute);


/****************************************************************************
 * Início do laço para cálculo RK4 Para um disparo ocorrido a partir do NMM *
 *                                                                          *
 ****************************************************************************/

    inclinacao = 1.0; // Qualquer valor positivo só para iniciar. Valor será recalculado dentro do laço. Essa variável será comparada com o valor real medido ϕ durante a perícia de local de crime.
    inclinacao_lateral = 1.0; // Qualquer valor positivo para entrar no laço.
    
        //Segundo teste do while: se o disparo for descendente, o laço só para quando o projétil cair abaixo da altura medida, caso contrário, para qunaodo o projétil estiver acima da altura medida

        
    while ( (phi<0) ? ((inclinacao >= 0) || (projetil_1.y>altura)) : (projetil_1.y<altura) ){
        t1 = t + H;
//funcao pra fazer tudo - call by reference
//runge_kutta só recebe como argumento uma função f e um passo; <- na verdade tem que receber todos os argumentos
        projetil_1.x = projetil.x + runge_kutta(&pos_x, &projetil, &w, H, 0, 0)*H;
        projetil_1.y = projetil.y + runge_kutta(&pos_y, &projetil, &w, H, 0, 0)*H;
        projetil_1.z = projetil.z + runge_kutta(&pos_z, &projetil, &w, H, 0, 0)*H;

        inclinacao_lateral = atan (projetil.vz/projetil.vx);

        kappa = kappaSdensidade*densidade_ar(projetil.y);

        //Lembrar que, a cada passo, o azimute atual muda, pois muda a inclinacao_lateral.

        projetil_1.vx = projetil.vx + runge_kutta(&w_vx, &projetil, &w, H, kappa, g)*H;
        projetil_1.vy = projetil.vy + runge_kutta(&w_vy, &projetil, &w, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + runge_kutta(&w_vz, &projetil, &w, H, kappa, g)*H;


        inclinacao = atan (projetil_1.vy/projetil_1.vx);
        //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição do projétil.
        
        //Atualizações de variáveis.
        t = t1;
        projetil = projetil_1;

        //Se com arrasto, o projétil caiu antes de atingir a "altura" é porque para o dado theta, ele não subirá muito. Então deve-se incrementar "theta".
        if (projetil.y<0 ) {
        //    printf("\n%lf\t%lf",theta,projetil.x);
            theta = theta + PARADA;

            if (theta>0.785) { //Se theta estiver a 45 graus e ainda sim o projétil não atingir a altura, não há solução para o problema com os parâmetros fornecidos. O programa encerra-se.
                printf ("\nCom os parâmetros fornecidos, o projétil não atingiria a altura de impactação.");
                fclose(arquivo);
                goto FATALERROR;
            }
        goto A;
        }

    }
    n++;

    inclinacao_lateral = atan (projetil_1.vz/projetil_1.vx);
    //inclinacao lateral final (como não é critério de parada para o laço, só precisa ser calculado ao final do processo).
    t_total = t; //Guarda o tempo total para colocar na ordem inversa na gravação do arquivo.


/************************************
 * Etapa de correção do θ inicial   *
 * de atrito inicial.               *
 *                                  *
 ************************************/

    delta_phi = inclinacao - phi;

    if (fabs (delta_phi) >  PARADA){
        theta = ajustar_theta (inclinacao,phi,theta);
        goto A;
    } 

    
/********************************************************************************************************
 * Etapa de correção do AZ0 (Azimute Inicial)                                                           *
 *                                                                                                      *
 * Considerando que γ é a inclinação lateral medida (em relação ao Norte em sentido horário),           *
 * tal qual o φ para a angulação vertical, para descobrir o azimute inicial do disparo, é               *
 * necessário que a inclinaçao lateral final + o azimute - γ seja menor que um delta arbitrariamente    *
 * escolhido. Ou seja, como γ é fixo e a inclinação lateral varia com a simulação, vamos corrigir o     *
 * azimute (incrementando ou subtraindo) para obter um resultado menor que um delta escolhido.          *
 *                                                                                                      *
 ********************************************************************************************************/
    
    delta_inclinacao_lateral = projetil.azimute + inclinacao_lateral - gamma;
    //delta_inclinacao_lateral e gamma (γ) precisam ser declarados. gamma (γ) é uma variavel medida, como φ.

    if( fabs (delta_inclinacao_lateral) > PARADA ){
        projetil.azimute = ajuste_AZ(projetil.azimute, inclinacao_lateral, gamma);
        goto B;
    }


#if DEBUG   /*Cálculo de duração de funções*/
printf("\n\n\nTEMPO GASTO NO PRIMEIRO LAÇO DE CALCULOS:\nt = %f segundos\n\n\n",  ((double) (clock() - ciclos_cpu))/CLOCKS_PER_SEC);
#endif
    
/****************************************************************************************************
 * A primeira simulação terminada forneceu as distâncias máximas possíveis para origem do disparo,  *
 * pois foi considerado que o disparo ocorreu ao nível médio do mar NMM (projetil.y "inicial" = 0)  *
 * Com isso, tem-se um novo limitante para futuras iterações (distância máxima percorrida).         *
 * Essa distância máxima percorrida é aproximadamente o donwrange máximo alcançado pelo projétil.   *
 ****************************************************************************************************/

    origemNMM[0] = (180/M_PI)*latitude - distLatGraus (&projetil);  //latitude ao Nível do Mar
    origemNMM[1] = longitude - distLongGraus (&projetil,latitude);  //longitude ao Nível do Mar
    double downrangeMax = projetil.x; // pra subtrair do Downrange do NMM até o prédio e ficar somente a distância entre o prédio e a impactacao.
    latitude_disparo = origemNMM[0];
    longitude_disparo = origemNMM[1];
    //Caso haja uma edificação, esse valor será atualizado se o projétil partir dessa edificação.

/************************************
 * Fim das primeiras simulações     *
 * considerandoo projétil partindo  *
 * do solo.                         *
 ************************************/

    /* Questionamento se havia edificações no caminho do projétil */
    printf("Existe alguma edificacao entre o impacto e a possível origem do disparo no solo de onde possa ter partido o tiro?\n");
    printf("\t\tLatitude\tLongitude\nImpacto\t\t%f,\t%f\nOrigem (NMM)\t%f,\t%f\n",(180/M_PI)*latitude,longitude,origemNMM[0], origemNMM[1]);
    printf("\n1 - SIM. 2 - NAO: ");

    if (getchar() == '1'){ /* Em ASCII: 49 */
            haviaEdf = 1;
#if DEBUG //, Holiday -8.123596, -34.898176 60m   , 
latitudeEdf=-8.123596;
longitudeEdf= -34.898176;
alturaEdf=60;
distanciaPredio_Impacataco = haversine (180*latitude/M_PI,longitude,latitudeEdf, longitudeEdf);

#else
            printf("\nInsira as coordenadas da edificação mais próximas da região da impactação.");
            printf("\nLatitude da edificação: ");
            ("%lf", &latitudeEdf);
            printf("Longitude da edificação: ");
            scanf("%lf", &longitudeEdf);
            printf("Altura da edificação: ");
            scanf("%lf", &alturaEdf);
            distanciaPredio_Impacataco = haversine (180*latitude/M_PI,longitude,latitudeEdf, longitudeEdf);
    #endif
    } else {
        printf("\nComo não há outras edificações no caminho:\n");
        altura_disparo=0;
        goto ULTIMA;
    }

#if DEBUG   /*Cálculo de duração de funções*/
ciclos_cpu = clock();
#endif

    C: //Ponto de partida do goto após correção do ângulo θ = theta/Phi e Altura para o sistema com edificação.

    //Condições iniciais para variáveis com edificação no caminho:
    
    t=0.0;
    projetil.x=0.0;
    projetil.y=altura_disparo;
    projetil.z=0.0;
    projetil.vx=v*cos(theta);
    projetil.vy=v*sin(theta);
    projetil.vz=0.0;   //Downrange continua sendo no eixo x;
    //projetil.azimute será atualizado nos laços.
    projetil_1 = projetil;
    
    //Recálculo do vento pela variação do azimute inicial.
    //Mudança de base de Norte-Leste para X e Z
    w.x = vN*cos(projetil.azimute) + vE*sin(projetil.azimute);
    w.y = 0.0;
    w.z = -vN*sin(projetil.azimute) + vE*cos(projetil.azimute);

/********************************************************************************************************************
 * Início do laço para cálculo RK4 para um disparo ocorrido de uma edificação mais próxima da região de impactação  *
 *                                                                                                                  *
 ********************************************************************************************************************/

    while ( projetil_1.x < downrangeMax){
        //A única coisa que importa após o edf. é a distância percorrida ser menor que a distância entre a edificação e a impactação.
        
        t1 = t + H;
        projetil_1.x = projetil.x + runge_kutta(&pos_x, &projetil, &w, H, 0, 0)*H;
        projetil_1.y = projetil.y + runge_kutta(&pos_y, &projetil, &w, H, 0, 0)*H;
        projetil_1.z = projetil.z + runge_kutta(&pos_z, &projetil, &w, H, 0, 0)*H;

        inclinacao_lateral = atan (projetil.vz/projetil.vx);

        kappa = kappaSdensidade*densidade_ar(projetil.y);

        //Lembrar que, a cada passo, o azimute atual muda, pois muda a inclinacao_lateral.

        projetil_1.vx = projetil.vx + runge_kutta(&w_vx, &projetil, &w, H, kappa, g)*H;
        projetil_1.vy = projetil.vy + runge_kutta(&w_vy, &projetil, &w, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + runge_kutta(&w_vz, &projetil, &w, H, kappa, g)*H;
        inclinacao = atan (projetil_1.vy/projetil_1.vx);
        //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição do projétil.

        if (haviaEdf) {

            if (projetil.x > (downrangeMax-distanciaPredio_Impacataco) ) { 
            //essa condicao aumenta o erro se o usuário colocar uma edificao fora da trajetória
            //Essa condicao testa para ver se o projétil já passou da edificação. Nesse caso, testa (abaixo) para saber se o projétil passou por cima da edificação.

                if(projetil.y < alturaEdf){

                    t=0.0;
                    downrangeMax = downrangeMax-projetil.x; // subtraindo a posição atual do projétil, temos a distância da edificação até a impactacao.
                    projetil.x=0.0; // Eixo "x" é o Downrange.
                    /*Gravação da variável altura_disparo.*/
                    altura_disparo=projetil.y;
                    projetil.y = altura_disparo; //Nada muda em projetil.y, ou seja, continua com a mesma altura. 
                    projetil.z=0.0;
                    /* Atualização do theta a partir da edificação.*/
                    theta = inclinacao;
                    projetil.vx=v*cos(theta);
                    projetil.vy=v*sin(theta);
                    projetil.vz=0.0; // Seria possível já iniciar projetil.vz com a inclinacao_lateral, mas deixa recalcular... menos trabalho.
                    projetil_1=projetil;
                    latitude_disparo = latitudeEdf;
                    longitude_disparo = longitudeEdf;
                    haviaEdf = 0; // == 0 Interrompe os testes para esses cálculos.

                } else {
                    printf ("O projétil passou por cima da edificação.");
                    altura_disparo=0;
                    latitude_disparo = origemNMM[0];
                    longitude_disparo = origemNMM[1];
                    goto ULTIMA;
                }
            }
        }

        //Atualizações de variáveis.
        t = t1;
        projetil = projetil_1;

        }
    n++;
    inclinacao_lateral = atan (projetil_1.vz/projetil_1.vx);
    //inclinacao lateral final (como não é critério de parada para o laço, só precisa ser calculado ao final do processo).
    t_total = t; //Guarda o tempo total para colocar na ordem inversa na gravação do arquivo.


/************************************
 * Etapa de correção do θ inicial   *
 * de atrito inicial.               *
 *                                  *
 ************************************/

    delta_phi = inclinacao - phi;
    
    if (fabs (delta_phi) >  PARADA){
        theta = ajustar_theta (inclinacao,phi,theta);
        goto A;
    } 

/****************************************************************
 * Etapa de correção da altura inicial do disparo a partir      *
 * do edificio.                                                 *
 *                                                              *
 * IMPORTANTE: Após a correção realizada no theta, o programa   *
 * ajustará a altura de disparo, mas note que é possível que:   *
 * com v0 mais alto, o projétil precise partir de uma altura    *
 * maior que a do próprio prédio.                               *
 * Quando foi inserido o prédio entre o solo e o disparo, o     *
 * projétil passava por esta região da edificacao com uma       *
 * velocidade MENOR do que V0, com isso conseguiria curvar e    *
 * atingir na região de impactação com a inclinação medida.     *
 ****************************************************************/

    delta_y = fabs (projetil.y - altura);

    if ( delta_y > 0.01){
        if (projetil.y > altura) {
            altura_disparo -= delta_y/2 + 0.001; //em delta_y já foi aplicado fabs.
            goto C;
        } else {
            altura_disparo += delta_y/2 + 0.001;
            goto C;
        }
    }

    /*Checagem para ver se o projétil parte acima da edificação*/
    if (altura_disparo > alturaEdf + 1.5) { /*+1.5 porque uma pessoa, em tese, pode subir no topo e disparar. Cuidado porque o tiro pode vir de antes do prédio. Analisar com calma.*/
        printf("\nATENÇÃO!\n");
        printf("\nPara que o projétil atinja uma altura de %.2lf m com inclinação de %.0lfº, disparado a %.0lf m de distância da impactação (Coordenadas: %lf ,%lf), precisaria ter sido disparado a uma altura de %.2lf m do solo (considerando que sua velocidade inicial é %.0lf m/s), ou seja, mais alto que a suposta edificação de origem, logo, não pode ter partido desta edificação.",altura,180*phi/M_PI,distanciaPredio_Impacataco,latitudeEdf,longitudeEdf,altura_disparo,v);
        printf("\n\nReveja as condições inciais e reinicie o programa.\n\n");
        fclose(arquivo);
        goto FATALERROR;
    }

/********************************************************************************************************
 * Etapa de correção do AZ0 (Azimute Inicial)                                                           *
 *                                                                                                      *
 * Considerando que γ é a inclinação lateral medida (em relação ao Norte em sentido horário),           *
 * tal qual o φ para a angulação vertical, para descobrir o azimute inicial do disparo, é               *
 * necessário que a inclinaçao lateral final + o azimute - γ seja menor que um delta arbitrariamente    *
 * escolhido. Ou seja, como γ é fixo e a inclinação lateral varia com a simulação, vamos corrigir o     *
 * azimute (incrementando ou subtraindo) para obter um resultado menor que um delta escolhido.          *
 *                                                                                                      *
 ********************************************************************************************************/
    
    delta_inclinacao_lateral = projetil.azimute + inclinacao_lateral - gamma;
    //delta_inclinacao_lateral e gamma (γ) precisam ser declarados. gamma (γ) é uma variavel medida, como φ.
    if( fabs (delta_inclinacao_lateral) > PARADA ){
        projetil.azimute = ajuste_AZ(projetil.azimute, inclinacao_lateral, gamma);
        goto B;
    }


#if DEBUG   /*Cálculo de duração de funções*/
printf("\n\n\nTEMPO GASTO NO SEGUNDO LAÇO DE CALCULOS:\nt = %f segundos\n\n\n",  ((double) (clock() - ciclos_cpu))/CLOCKS_PER_SEC);
#endif    


/************************************************************************
 * Recalculo a partir da última correção apenas para gravar no arquivo  *
 *                                                                      *
 ************************************************************************/    

    ULTIMA: //Ponto de partida do goto quando não há edificações no percurso do projétil.
    

    t=0.0;
    projetil.x=0.0; // Eixo "x" é o Downrange.
    projetil.y=altura_disparo;
    projetil.z=0.0;
    projetil.vx=v*cos(theta);
    projetil.vy=v*sin(theta);
    projetil.vz=0.0;
    projetil_1=projetil;
    
    inclinacao = 1.0;
    // Valor arbitrário positivo para poder entrar no critério do laço. Valor será recalculado dentro do laço.

/************************************
 * Impressão das condições iniciais *
 *                                  *
 ************************************/

    fprintf(arquivo,"A posição e altura do projétil dadas a partir do momento do disparo na seguinte ordem:"
                    "\nTempo\tLatitude\t\t\tLongitude\t\tAltura:\n");

    /* O último laço para gravar os resultados precisa apenas obedecer o downrangeMax com as devidas condições iniciais. */
    while ( projetil_1.x < downrangeMax){
        t1 = t + H;
        projetil_1.x = projetil.x + runge_kutta(&pos_x, &projetil, &w, H, 0, 0)*H;
        projetil_1.y = projetil.y + runge_kutta(&pos_y, &projetil, &w, H, 0, 0)*H;
        projetil_1.z = projetil.z + runge_kutta(&pos_z, &projetil, &w, H, 0, 0)*H;

        inclinacao_lateral = atan (projetil.vz/projetil.vx);
        
        kappa = kappaSdensidade*densidade_ar(projetil.y);
 
/********************************
 * Imprimir alguma informação   *
 * relevante para análise.      *
 *                              *
 ********************************/
#if DEBUG
fprintf(debug,"%lf,%lf\tAltura: %lf m\n",latitudeEdf + distLatGraus (&projetil_1), longitudeEdf + distLongGraus(&projetil_1,latitude),projetil.y);
#endif

        projetil_1.vx = projetil.vx + runge_kutta(&w_vx, &projetil, &w, H, kappa, g)*H;
        projetil_1.vy = projetil.vy + runge_kutta(&w_vy, &projetil, &w, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + runge_kutta(&w_vz, &projetil, &w, H, kappa, g)*H;
        inclinacao = atan (projetil_1.vy/projetil_1.vx);

        //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição específica do projétil.

        t = t1;
        projetil = projetil_1;        
       
        fprintf(arquivo, "%.3lf\t%.12lf, %.12lf\t %lf m\n", t1, latitude_disparo + distLatGraus (&projetil_1), longitude_disparo + distLongGraus(&projetil_1,latitude),projetil_1.y);
        
    }

    if (projetil_1.azimute > 2*M_PI) projetil_1.azimute = projetil_1.azimute - 2*M_PI; //Para deixar entre 0º e 360º, porque ficou maior que 360º.
    if (projetil_1.azimute < 0) projetil_1.azimute = projetil_1.azimute + 2*M_PI;      //Para deixar entre 0º e 360º, porque ficou menor que 0º.
    
/************************************
 * Considerações Finais e           *
 * exibição em tela do resultado    *
 ************************************/    

    printf("\nConsiderando um sistema com arrasto, os cálculos terminaram com os seguintes valores:"
           "\nForam efetuados %d cálculos de trajetória."
           "\nDownrange Total = %.3lf m."
           "\nAltura de impactação = %.3lf m."
           "\nDesvios para %s = %.3lf m."
           "\nÂngulo θ (inicial) do disparo = %.2lfº."
           "\nÂngulo ϕ (ao final da simulação) = %.2lfº."
           "\nAzimute inicial do disparo = %.2lfº.\n",n,projetil.x,projetil.y,(projetil.z<0 ? "esquerda" : "direita"), fabs(projetil.z), 180*theta/M_PI,180*inclinacao/M_PI, 180*projetil.azimute/M_PI);
    printf("\nA trajetória teve outro desvio devido ao spindrift de, aproximadamente, %.0f cm para %s, não incluidos nos cálculos.\n", spindrift(t), dextrogiro ? "direta" : "esquerda");
    
    printf("\nO projétil partiu, aproximadamente, das coordenadas: %lf N/S, %lf L/O, a uma altura de %.2lf m.\n",latitude_disparo, longitude_disparo, altura_disparo); //latitude_disparo e longitude_disparo foram calculados ao fim do segundo laço.
 
    velocidade_final = sqrt (pow(projetil.vx,2)+pow(projetil.vy,2)+pow(projetil.vz,2)); //Módulo nas três componentes.
    printf("\nO projétil levou cerca de %.1f segundos do momento do disparo à impactação.\nPossuía velocidade final de %.2f m/s e energia cinética de %.2f J.\n",t,velocidade_final,0.5*m*pow(velocidade_final,2));
    
/************************************************
 * Dados para comparação da energia cinética    *
 *                                              *
 ************************************************/
    printf("\nEnergia cinética de alguns projéteis para comparação:\n");
    printf("Calibre \tEnergia (J)\n");
    printf(".25 AUTO\t87     \n");
    printf(".32 AUTO\t175    \n");
    printf(".380 AUTO EXPO  259\n");
    printf(".38 SPL CHOG    271\n");
    printf("9x19mm (124gr)  459\n");
    printf(".40 EXPO Gold   568\n");
    printf(".357 Mag        724\n");
    printf(".454 Casull     2531\n");

    
    fclose(arquivo);
    
#if DEBUG
    fclose(debug);
#endif
    
    FATALERROR:
    return 0;
}
