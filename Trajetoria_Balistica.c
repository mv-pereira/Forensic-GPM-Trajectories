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


#define OMEGA 0.000072921   //Taxa de rotação da terra em "rad/s".
#define PARADA 0.0174533    // Critério de parada para ajuste de angulação. 0.0174533 rad = 1º.
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

/****************************************************************************
 *                        Spindrift aproximado                              *
 * Valores estimados para projéteis subsônicos a partir de                  *
 * https://theoverwatch.wixsite.com/theoverwatch/post/spin-drift.           *
 * sg = 2.42563 para projétil subsônico (300 Blackout using subsonic ammo)  *
 * sg = 1.52551 para projétil supersônico .308                              *
 ****************************************************************************/
double spindrift (double tempo){
    double sg,drift;
    sg = 2.42563;
    drift = 1.25*(sg+1.2)*pow(tempo,1.83);  //Resultado em polegada.
    return 2.54*drift;                      //Resultado em centímetro.
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

double kvx (struct prjt *projetil, struct vento *w, double inclinacao_lateral, double h, double kappa){
    double k,kvx1,kvx2,kvx3,kvx4;
    projetil->azimute += inclinacao_lateral;
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


double kvy (struct prjt *projetil, struct vento *w, double inclinacao_lateral, double h, double kappa, double g){
    double k,kvy1,kvy2,kvy3,kvy4;
    projetil->azimute += inclinacao_lateral;
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
 *                                            projetil_1.vz1.     *
 *                                                                *
 ******************************************************************/

double w_vz (double k, struct prjt *projetil, double correcao, struct vento *w){
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vz +correcao -w->z) + 2*OMEGA*( (projetil->vx +correcao -w->x)*(sin (projetil->latitude)) -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(cos (projetil->azimute)) )   );
}


double kvz (struct prjt *projetil, struct vento *w, double inclinacao_lateral, double h, double kappa){
    double k,kvz1,kvz2,kvz3,kvz4;
    projetil->azimute += inclinacao_lateral;
    kvz1 = w_vz(kappa, projetil, 0, w);
    kvz2 = w_vz(kappa, projetil, kvz1*(h/2), w);
    kvz3 = w_vz(kappa, projetil, kvz2*(h/2), w);
    kvz4 = w_vz(kappa, projetil, kvz3*(h/2), w);
    k = (1/6.0)*(kvz1 + 2*kvz2 + 2*kvz3 + kvz4);
    return k;
}

/****************************************************
 * Função Principal: Trajetória Balística           *
 *                                                  *
 ****************************************************/

int main(){
    setlocale(LC_ALL, "Portuguese"); //Utilizando caracteres e acentuação da língua portuguesa.

    FILE *arquivo;
    arquivo = fopen("data","w");
    
#if DEBUG //DEBUG para geração de arquivos com informações pertinentes.
        FILE *debug;
        debug = fopen("debug","w"); 
#endif

    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        exit(1);
    }
    
    int n;                              //Contador.
    int descendente = 1, destrogiro;    //variável "booleana" para critério do laço, e se o projétil é destrogiro ou levogiro.
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
    
    struct prjt projetil, projetil_1; //projetil_1 é o projétil para n+1 | Definição dessa struct no começo do programa.
    struct vento w;                     //Definição da struct do vento.

#if DEBUG //Valores padrão apenas para comparação

altura = 89;
descendente = 0;
phi = 16*M_PI/180;
gamma = 73*M_PI/180;
v = 271;
m=(10.24)/1000.0;
a=M_PI*powl((8.82/1000.0),2)/4; // A divisão por quatro leva em conta o raio.
destrogiro=1;
c = 0.23;
vento = 7.56/3.6;
ang = (200+180.0)*M_PI/180.0;
latitude = -8.108258*M_PI/180.0;
longitude = -34.892634;
azimute = 73*M_PI/180.0;
g = 9.780327*(1+0.0053024*sin(latitude)*sin(latitude) - 0.0000058*sin(2*latitude)*sin(2*latitude)); //Açeleração da gravidade na latitude. (em m/s^2)
printf ("\n*\t*\tDEBUG Ativado.\t*\t*\n\t\tValores prefixados.\n\nPara sair da função DEBUG, mudar a definição de DEBUG para 0 no cabeçalho do programa e recompilar.\n\n");

#else

/********************************
 * Características do projetil  *
 ********************************/
    
    printf("Digite a altura (em metros) da impactação em relação ao disparo: ");
    scanf("%lf", &altura);
    
    printf("A impactação ocorreu em trajetória descendente?\n1 - SIM\t2 - NÃO.\n");
    scanf("%d", &descendente);
    (descendente != 1) ? descendente = 0 : 1;
    //Forçar descdente ser 1 (True) ou 0 (False).

    printf("\nDigite o ângulo ϕ com a horizontal (em °): ");
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
    
    printf("\nO projetil é destrogiro ou levogiro? 1 - Destrogiro.\t2 - Levogiro.\n");
    scanf("%d", &destrogiro);
    (destrogiro != 1) ? destrogiro = 0 : 1;
 
    printf("\nDigite o módulo da velocidade do vento (em km/h): ");
    scanf("%lf", &vento);
    vento = vento/3.6; // Velocidade do vento precisa ser em m/s.

    printf("\nDigite a direção do vento em relação ao Norte (em °): ");
    scanf("%lf", &ang);
    ang = (ang+180.0)*M_PI/180.0; // A direção dos ventos é dada de onde ele vem, não para onde sopra.
    
    printf("\nDigite a latitude decimal do disparo (em °), Valor precisa ser negativo, se o disparo ocorrer Hemisfério Sul: ");
    scanf("%lf", &latitude);
    latitude = latitude*M_PI/180.0; //grau para radianos

    printf("\nDigite a longitude decimal do disparo (em °): ");
    scanf("%lf", &longitude);
    
    printf("\nDigite o azimute estimado do disparo (em °), medido em sentido horário em relação ao Norte: ");
    scanf("%lf", &azimute);
    azimute = azimute*M_PI/180.0; //grau para radianos


#endif    

/****************************************************************************
 * Velocidade do vento é dada a partir de onde ele sopra e será decomposto  *
 * primeiramente nos eixos Norte e Leste e, posteriormente, nos eixos x     *
 * (Downrange) e z (desvio lateral)                                         *
 ****************************************************************************/

    vN = vento*cos(ang);
    vE = vento*sin(ang);

    /* TIRAR ESSE PEDAÇO
    w.x = vN*cos(azimute) + vE*sin(azimute); //mudança de base de Norte-Leste para X e Z
    w.y = 0.0;
    w.z = -vN*sin(azimute) + vE*cos(azimute); //mudança de base de Norte-Leste para X e Z    
    */

    g = 9.780327*(1+0.0053024*pow(sin(latitude),2) - 0.0000058*pow (sin(2*latitude),2)); //Açeleração da gravidade na latitude. (em m/s^2) nvl do mar

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
    
    kappaSdensidade = c*a/(2.0*m); //Constante de arrasto/ro. Note que a massa já entra na constante. Densidade do ar será calculado nas funcoes de vx,vy,vz...
    projetil.latitude=latitude;
    projetil.azimute =azimute;
    n = 0;                          //Contador para registrar quantas vezes toda a trajetória foi calculada.
    
/********************************************************
 * Condições Iniciais: Dados dispostos considerando um  *
 *                     disparo a um metro de altura     *
 *                     a 0 grau (projetil.vy=0).        *
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

    
/************************************
 * Início do laço para cálculo RK4  *
 *                                  *
 ************************************/

    inclinacao = 1.0; // Qualquer valor positivo só para iniciar. Valor será recalculado dentro do laço. Essa variável será comparada com o valor real medido ϕ durante a perícia de local de crime.
    inclinacao_lateral = 1.0; // Qualquer valor positivo para entrar no laço.
    
    //((gamma > 0) ? (projetil_1.y>altura) : (projetil_1.y<altura))
        //Segundo teste do while: se o disparo for descendente, o laço só para quando o projétil cair abaixo da altura medida, caso contrário, para qunaodo o projétil estiver acima da altura medida
    
    while ( descendente ? ((inclinacao > 0) || (projetil_1.y>altura)) : (projetil_1.y<altura)){
        t1 = t + H;
        projetil_1.x = projetil.x + pos(projetil.vx,H)*H;
        projetil_1.y = projetil.y + pos(projetil.vy,H)*H;
        projetil_1.z = projetil.z + pos(projetil.vz,H)*H;

        inclinacao_lateral = atan (projetil.vz/projetil.vx);
        
        kappa = kappaSdensidade*densidade_ar(projetil.y);

        projetil_1.vx = projetil.vx + kvx(&projetil, &w, inclinacao_lateral, H, kappa)*H;   //Lembrar que, a cada passo, o azimute atual muda, pois muda a inclinacao_lateral.
        projetil_1.vy = projetil.vy + kvy(&projetil, &w, inclinacao_lateral, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + kvz(&projetil, &w, inclinacao_lateral, H, kappa)*H;
        inclinacao = atan (projetil_1.vy/projetil_1.vx);
        //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição do projétil.
        
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
    delta_phi = fabs (inclinacao) - phi;
    
    if ( fabs (delta_phi) > PARADA/10){ // Critério de parada. Se o Δϕ for maior que 0.1º (0.00174533 rad), soma ou subtrai 0.01º.
    

    if ( delta_phi > 0 ){               //Projétil terminou com angulação maior que o φ medido, indicando que o disparo foi mais baixo.
        theta=theta-PARADA/100;         //Diminue θ inicial.
        goto A;
        }
    else{                               //Projétil terminou com angulação menos que o φ medido, indicando que o disparo foi mais alto.
        theta=theta+PARADA/100;         //Aumenta θ inicial.
        goto A;
        }
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
    if ( fabs (delta_inclinacao_lateral) > PARADA/10){      // Critério de parada. Se for maior que 0.1º (0.00174533 rad), soma ou subtrai 0.01º.

    if (delta_inclinacao_lateral > 0){                      //Projétil terminou com azimute maior que o γ medido, devendo reduzir azimute inicial.
        projetil.azimute = projetil.azimute - PARADA/100;   //Diminue Azimute inicial.
        goto B;
        }
    else{                                                   //Projétil terminou com azimute maior que o γ medido, devendo reduzir azimute inicial.
        projetil.azimute = projetil.azimute + PARADA/100;   //Diminue Azimute inicial.
        goto B;
        }
    }
    
/************************************
 * Recalculo a partir da última     *
 * correção apenas para gravar      *
 * no arquivo                       *
 ************************************/    

    t=0.0;
    projetil.x=0.0; // Eixo "x" é o Downrange.
    projetil.y=0.0;
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

    fprintf(arquivo,"Os dados serão gravados (distâncias em m e velocidades em m/s)"
                    " na ordem:\ntempo    L/O      Altura    N/S      V(L/O)      V(y)      V(N/S)   Downrange:\n");
    fprintf(arquivo,"%lf %lf %lf %lf %lf %lf %lf %lf\n",t,projetil.x,projetil.y,projetil.z,projetil.vx,projetil.vy,projetil.vz, projetil.x);

    while ( descendente ? ( (inclinacao > 0) || (projetil_1.y>altura) ) : (projetil_1.y<altura)){
        t1 = t + H;
        projetil_1.x = projetil.x + pos(projetil.vx,H)*H;
        projetil_1.y = projetil.y + pos(projetil.vy,H)*H;
        projetil_1.z = projetil.z + pos(projetil.vz,H)*H;

        inclinacao_lateral = atan (projetil.vz/projetil.vx);
        
        kappa = kappaSdensidade*densidade_ar(projetil.y);
 

#if DEBUG
/********************************
 * Imprimir alguma informação   *
 * relevante para análise.      *
 *                              *
 ********************************/
fprintf(debug,"%f\n",projetil.vx);
#endif

        projetil_1.vx = projetil.vx + kvx(&projetil, &w, inclinacao_lateral, H, kappa)*H;   //Lembrar que, a cada passo, o azimute atual muda, pois muda a inclinacao_lateral.
        projetil_1.vy = projetil.vy + kvy(&projetil, &w, inclinacao_lateral, H, kappa, g)*H;
        projetil_1.vz = projetil.vz + kvz(&projetil, &w, inclinacao_lateral, H, kappa)*H;
        inclinacao = atan (projetil_1.vy/projetil_1.vx);

        //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição específica do projétil.

        t = t1;
        projetil = projetil_1;        
       
        fprintf(arquivo,"%lf %lf %lf %lf %lf %lf %lf %lf\n",t1,(projetil_1.x*cos(projetil_1.azimute) - projetil_1.z*sin(projetil_1.azimute) ),projetil_1.y,(projetil_1.x*sin(projetil_1.azimute) + projetil_1.z*cos(projetil_1.azimute) ),projetil_1.vx,projetil_1.vy,projetil_1.vz,projetil_1.x); //x e z rotacionados para Norte-Leste
        
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
           //"\nÂngulo ϕ de impactação = %.2lfº."
           "\nÂngulo θ (inicial) do disparo = %.2lfº."
           "\nAzimute inicial do disparo = %.2lfº.\n",n,projetil.x,projetil.y,(projetil.z<0 ? "esquerda" : "direita"), fabs(projetil.z), 180*theta/M_PI, 180*projetil.azimute/M_PI);
    printf("\nA trajetória teve outro desvio devido ao spindrift de, aproximadamente, %.0f cm para %s, não incluidos nos cálculos.\n", spindrift(t), destrogiro ? "direta" : "esquerda");

/****************************************************************
 * Cálculo para imprimir na tela as coordenadas geográficas do  *
 * disparo considerando um tiro efetuado a partir do solo!      *
 * "distLatGraus" é a "distância" em graus entre a latitude do  *
 * disparo e a latitude da impactação.                          *
 * "distLongGraus" é a "distância" em graus entre a latitude do *
 * disparo e a latitude da impactação.                          *
 *                                                              *
 * Constante "0.000008983152098" tem unidades de º/m.           *
 ****************************************************************/

    double distLatGraus, distLongGraus;
    distLatGraus = 0.000008983152098*(projetil_1.x*cos(projetil_1.azimute) - projetil_1.z*sin(projetil_1.azimute));
    distLongGraus = 0.000008983152098*cos(latitude)*(projetil_1.x*sin(projetil_1.azimute) + projetil_1.z*cos(projetil_1.azimute) );
    
    printf("\nConsiderando um disparo a partir do solo, o projétil partiu, aproximadamente, das coordenadas: %f N/S, %f L/O.\n",(180/M_PI)*latitude- distLatGraus, longitude - distLongGraus);
    printf("\nO projétil levou cerca de %.1f segundos do momento do disparo à impactação.\nPossuía velocidade final de %.2f m/s e energia cinética de %.2f J.\n",t,projetil.vx,0.5*m*pow(projetil.vx,2));
    
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
    
    return 0;
}
