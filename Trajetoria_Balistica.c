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



#include <locale.h> //Utilizando caracteres e acentuação da língua portuguesa.
#include <stdlib.h> //Para função exit() na condicional da abertura do arquivo;
#include <stdbool.h>
#include <time.h>   //Para velocidade de funções - DEBUG
#include "projetil.h"


#define H 0.0001            //passo da iteração do Runge-Kutta.
#define DEBUG 0

double calcularg(double latitude);

void lerDadosDoArquivo(char *nomeArquivo, struct impactacao *impacto, struct disparo *tiro, struct prjt *projetil, struct vento *w, double *g);

void lerDadosEdificacao(char *nomeArquivo, struct edificacao *edificio);

void atualizarProjetil(struct prjt *projetil, struct vento *w, double g, double azimuteTiro);

bool correcaoTheta(struct prjt *projetil, struct disparo *tiro, const double impacto_phi);

bool correcaoAzimute(struct prjt *projetil, struct disparo *tiro, const double impacto_azimute);

bool correcaoAltura(struct prjt *projetil, struct disparo *tiro, const double impacto_altura);

void inicializarProjetilEW(double *t, struct prjt *projetil, const struct disparo *tiro, struct vento *w, const double yInicial);

void inicializarTiro(const struct prjt *projetil, const struct impactacao *impacto, struct disparo *tiro, const struct vento *w, const struct edificacao *edificio, enum origem_disparo origem);

double movimentoProjetil(int *n, struct prjt *projetil, struct impactacao *impacto, struct disparo *tiro, struct vento *w, struct edificacao *edificio, bool calcular_Edf);

int main(){
    setlocale(LC_ALL, "Portuguese"); //Utilizando caracteres e acentuação da língua portuguesa.

    // FILE *arquivo;
    // arquivo = fopen("data","w");

    // if (arquivo == NULL) {
    //     printf("Erro ao abrir o arquivo!\n");
    //     exit(1);
    // }

    struct prjt projetil;               //Estrutura do projétil.
    struct vento w;                     //Definição da struct do vento.
    struct edificacao edificio;         //Estrutura da edificação.
    struct impactacao impacto;          //Estrutura da Impactação.
    struct disparo tiro;                //Estrutura do Tiro.

    int n;                              //Contador.
    bool calcular_Edf = false;               //variável "booleana" para cálculo se disparo ocorreu em edificação.
    int dextrogiro;                     //variável "booleana" para critério do laço, e se o projétil é dextrogiro ou levogiro.
    double xs[2];                       //x+ e x- Raízes do sistema sem arrasto ("s" sem arrasto).
    double t,t_total,g;                 //Tempo, Tempo total e aceleração da gravidade.
    double downrangeMax;                //máxima distância alcançada pelo projétil.
    double distanciaPredio_Impacataco;  //Distancia entre a edificacao e a impactacao a ser calculada caso haja edificações.
    double velocidade_final;            //Velocidade do projétil na impactactação.

    bool corrigirTheta = true;          // Variáveis definidas para as correções de ângulação e altura do disparo
    bool corrigirAzimute = true;
    bool corrigirAltura = true;

    lerDadosDoArquivo("input.txt", &impacto, &tiro, &projetil, &w, &g);

/************************************
 * Início do cálculo SEM arrasto    *
 * para estimativa do θ inicial     *
 * a ser utilizado nas iterações    *
 ************************************/    

    tiro.theta = arcsec (sqrtl(pow(tiro.velocidade,2)*sec(impacto.phi)*sec(impacto.phi)/(pow(tiro.velocidade,2)-2*g*impacto.altura)));
    printf("\nO ângulo θ do início do disparo com a horizontal considerado a partir do solo e em um sistema sem arrasto ou outras correcoes vale: %.2lf°\n",tiro.theta*180/M_PI);

    xs[0] = (powl(cos(impacto.phi),2)*(pow(tiro.velocidade,2) - 2*g*impacto.altura)*( -fabs(tan(impacto.phi)) + sqrtl(pow(tiro.velocidade,2)*powl(sec(impacto.phi),2)/(pow(tiro.velocidade,2)-2*g*impacto.altura) - 1) ) ) /g ;
    xs[1] = (powl(cos(impacto.phi),2)*(pow(tiro.velocidade,2) - 2*g*impacto.altura)*( +fabs(tan(impacto.phi)) + sqrtl(pow(tiro.velocidade,2)*powl(sec(impacto.phi),2)/(pow(tiro.velocidade,2)-2*g*impacto.altura) - 1) ) ) /g ;

    printf("\nConsiderando um sistema sem arrasto: x- = %.3lf m e x+ = %.3lf m\n",xs[0],xs[1]);
    printf("Onde x- é a menor, e x+ a maior distância que o disparo foi efetuado a partir de um ângulo |θ|.\n");
    
/********************************
 * Fim do cálculo sem Arrasto   *
 * onde foi obtido um θ inicial *
 ********************************/

/****************************************************************************
 * Velocidade do vento é dada a partir de onde ele sopra e será decomposto  *
 * primeiramente nos eixos Norte e Leste e, posteriormente, nos eixos x     *
 * (Downrange) e z (desvio lateral)                                         *
 ****************************************************************************/

    w.norte = w.velocidade*cos(w.direcao);
    w.leste = w.velocidade*sin(w.direcao);

    n = 0;  //Contador para registrar quantas vezes toda a trajetória foi calculada.

/****************************************************************
 * Condições Iniciais:                                          *
 * Dados dispostos considerando um disparo a 0 metros de altura *
 * e projetil.vy=tiro.velocidade*sin(tiro.theta)                *
 *                                                              *
 ****************************************************************/    

    // Principal função
    movimentoProjetil(&n, &projetil, &impacto, &tiro, &w, &edificio, calcular_Edf);
    inicializarTiro(&projetil, &impacto, &tiro, &w, &edificio, Nivel_do_Mar);


/************************************
 * Fim das primeiras simulações     *
 * considerandoo projétil partindo  *
 * do solo.                         *
 ************************************/

    /* Questionamento se havia edificações no caminho do projétil */
    printf("Existe alguma edificacao entre o impacto e a possível origem do disparo no solo de onde possa ter partido o tiro?\n");
    printf("\t\tLatitude\tLongitude\nImpacto\t\t%f,\t%f\nOrigem (NMM)\t%f,\t%f\n",(180/M_PI)*impacto.latitude,impacto.longitude,tiro.latitude, tiro.longitude);
    printf("\n1 - SIM. 2 - NAO: "); 
    fflush(stdin);
    char predio_no_caminho;
    scanf ("%c",&predio_no_caminho);
    if (predio_no_caminho == '1'){ /* Em ASCII: 49 */
            calcular_Edf = true;
    } else {
        printf("\nComo não há outras edificações no caminho:\n");
        tiro.altura=0;
        tiro.origem = Nivel_do_Mar;

    }
    // Cálculo com edifício
    t = movimentoProjetil(&n, &projetil, &impacto, &tiro, &w, &edificio, calcular_Edf);

/************************************
 * Considerações Finais e           *
 * exibição em tela do resultado    *
 ************************************/    

    printf("\nConsiderando um sistema com arrasto, os cálculos terminaram com os seguintes valores:"
           "\nForam efetuados %d iterações para cálculos de trajetória."
           "\nDownrange Total = %.3lf m."
           "\nAltura de impactação = %.3lf m."
           "\nDesvios para %s = %.3lf m."
           "\nÂngulo θ (inicial) do disparo = %.2lfº."
           "\nÂngulo ϕ (ao final da simulação) = %.2lfº."
           "\nAzimute inicial do disparo = %.2lfº.\n",n, projetil.x, projetil.y,(projetil.z<0 ? "esquerda" : "direita"), fabs(projetil.z), 180*tiro.theta/M_PI, 180*projetil.taxa_de_subida/M_PI, 180*tiro.azimute/M_PI);
    printf("\nA trajetória teve outro desvio devido ao spindrift de, aproximadamente, %.0f cm para %s, não incluidos nos cálculos.\n", spindrift(t), projetil.propriedades.rotacao == Dextrogiro ? "direta" : "esquerda");
    
    printf("\nO projétil partiu, aproximadamente, das coordenadas: %lf N/S, %lf L/O, a uma altura de %.2lf m, partindo %s.\n",tiro.latitude, tiro.longitude, tiro.altura, (tiro.origem == Nivel_do_Mar ? "ao nível do mar" : "de uma edificação" ) ); //latitude_disparo e longitude_disparo foram calculados ao fim do segundo laço.
 
    velocidade_final = sqrt (pow(projetil.vx,2)+pow(projetil.vy,2)+pow(projetil.vz,2)); //Módulo nas três componentes.
    printf("\nO projétil levou cerca de %.1f segundos do momento do disparo à impactação.\nPossuía velocidade final de %.2f m/s e energia cinética de %.2f J.\n",t,velocidade_final,0.5*projetil.propriedades.massa*pow(velocidade_final,2));
    
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

    return 0;
}





void lerDadosDoArquivo(char *nomeArquivo, struct impactacao *impacto, struct disparo *tiro, struct prjt *projetil, struct vento *w, double *g) {
    
    FILE *file = fopen(nomeArquivo, "r");
    char buffer[256];

    if (file == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        exit(1);
    }

    // Ignorando comentários e linhas em branco no início do arquivo
    while (fgets(buffer, sizeof(buffer), file)) {
        if (buffer[0] != '#' && buffer[0] != '\n') {
            break;
        }
    }

    sscanf(buffer, "%*[^:]: %lf", &impacto->altura);
    printf("Altura de Impacto: %lf m\n", impacto->altura);

    fscanf(file, "%*[^:]: %lf", &impacto->phi);
    printf("Ângulo Phi: %lf graus\n", impacto->phi);
    impacto->phi = impacto->phi*M_PI/180;

    fscanf(file, "%*[^:]: %lf", &impacto->azimute);
    printf("Ângulo Azimute: %lf graus\n", impacto->azimute);
    impacto->azimute = impacto->azimute*M_PI/180;

    fscanf(file, "%*[^:]: %lf", &impacto->latitude);
    printf("Latitude: %lf graus\n", impacto->latitude);
    impacto->latitude = impacto->latitude*M_PI/180.0;

    *g = calcularg(impacto->latitude); //Açeleração da gravidade na latitude. (em m/s^2)

    fscanf(file, "%*[^:]: %lf", &impacto->longitude);
    printf("Longitude: %lf\n", impacto->longitude);

    tiro->azimute = impacto->azimute;

    fscanf(file, "%*[^:]: %lf", &tiro->velocidade);
    printf("Velocidade Inicial: %lf m/s\n", tiro->velocidade);

    fscanf(file, "%*[^:]: %lf", &projetil->propriedades.massa);
    printf("Massa do Projétil: %lf g\n", projetil->propriedades.massa);
    projetil->propriedades.massa = projetil->propriedades.massa/1000.0;

    fscanf(file, "%*[^:]: %lf", &projetil->propriedades.diametro);
    printf("Diâmetro do Projétil: %lf mm\n", projetil->propriedades.diametro);
    projetil->propriedades.diametro = M_PI*powl((projetil->propriedades.diametro/1000.0),2)/4;  // A divisão por quatro leva em conta o raio.

    fscanf(file, "%*[^:]: %lf", &projetil->propriedades.coef_arrasto);
    printf("Coeficiente de Arrasto: %lf\n", projetil->propriedades.coef_arrasto);

    int rotacao;
    fscanf(file, "%*[^:]: %d", &rotacao);
    projetil->propriedades.rotacao = (rotacao == 1) ? Dextrogiro : Levogiro;
    printf("Rotação do Projétil: %d\n", rotacao);

    fscanf(file, "%*[^:]: %lf", &w->velocidade);
    w->velocidade = w->velocidade / 3.6; // Convertendo km/h para m/s
    printf("Velocidade do Vento: %lf m/s\n", w->velocidade);

    fscanf(file, "%*[^:]: %lf", &w->direcao);
    printf("Direção do Vento: %lf graus\n", w->direcao);
    w->direcao = (w->direcao+180.0)*M_PI/180.0;

    fclose(file);
}


void lerDadosEdificacao(char *nomeArquivo, struct edificacao *edificio) {
    FILE *file = fopen(nomeArquivo, "r");
    char buffer[256];

    if (file == NULL) {
        printf("Erro ao abrir o arquivo de edificação!\n");
        exit(1);
    }

    // Ignorando comentários e linhas em branco no início do arquivo
    while (fgets(buffer, sizeof(buffer), file)) {
        if (buffer[0] != '#' && buffer[0] != '\n') {
            break;
        }
    }

    sscanf(buffer, "%*[^:]: %lf", &edificio->latitude);
    printf("Latitude da Edificação: %lf graus\n", edificio->latitude);

    fscanf(file, "%*[^:]: %lf", &edificio->longitude);
    printf("Longitude da Edificação: %lf graus\n", edificio->longitude);

    fscanf(file, "%*[^:]: %lf", &edificio->altura);
    printf("Altura da Edificação: %lf metros\n", edificio->altura);

    fclose(file);
}


void atualizarProjetil(struct prjt *projetil, struct vento *w, double g, double azimuteTiro) {
    //Constante de arrasto/ro. Note que a massa já entra na constante. Densidade do ar será calculado nas funcoes de vx,vy,vz...
    //kappa - Variável da qual depende o arrasto. Há o cálculo inicial de tal kappa sem densidade e a cada altura (projetil.y) é calculada a densidade do ar.
    double kappaSdensidade = (projetil->propriedades.coef_arrasto/projetil->propriedades.massa)*(projetil->propriedades.diametro/2); //Ordem alterada de c*a/(2.0*m) para evitar Underflow.
    double kappa = kappaSdensidade * densidade_ar(projetil->y); // kappaSdensidade deve ser definido em algum lugar do seu código

    projetil->x += runge_kutta(&pos_x, projetil, w, H, kappa, g) * H;
    projetil->y += runge_kutta(&pos_y, projetil, w, H, kappa, g) * H;
    projetil->z += runge_kutta(&pos_z, projetil, w, H, kappa, g) * H;

    projetil->vx += runge_kutta(&w_vx, projetil, w, H, kappa, g) * H;
    projetil->vy += runge_kutta(&w_vy, projetil, w, H, kappa, g) * H;
    projetil->vz += runge_kutta(&w_vz, projetil, w, H, kappa, g) * H;

    projetil->taxa_de_subida = atan2(projetil->vy, projetil->vx);
    projetil->rumo = azimuteTiro + atan2(projetil->vz, projetil->vx);
}

bool correcaoTheta(struct prjt *projetil, struct disparo *tiro, const double impacto_phi) {
    double delta_phi = projetil->taxa_de_subida - impacto_phi;
    if (fabs(delta_phi) > PARADA) {
        tiro->theta = ajustar_theta(projetil->taxa_de_subida, impacto_phi, tiro->theta);
        return true; // Indica que ainda é necessária a correção
    }
    return false; // Não é mais necessário corrigir
}

bool correcaoAzimute(struct prjt *projetil, struct disparo *tiro, const double impacto_azimute) {
    double delta_rumo = projetil->rumo - impacto_azimute;
    if (fabs(delta_rumo) > PARADA) {
        tiro->azimute = ajuste_AZ(tiro->azimute, projetil->rumo, impacto_azimute);
        return true; // Indica que ainda é necessária a correção
    }
    return false; // Não é mais necessário corrigir
}

    /********************************************************************
     * Etapa de correção da altura inicial do disparo a partir          *
     * do edificio.                                                     *
     *                                                                  *
     * IMPORTANTE: Após a correção realizada no theta, o programa       *
     * ajustará a altura de disparo, mas note que é possível que:       *
     * com velocidade inicial do tiro mais alta, o projétil precise     *
     * partir de uma altura maior que a do próprio prédio.              *
     * Quando foi inserido o prédio entre o solo e o disparo, o         *
     * projétil passava por esta região da edificacao com uma           *
     * velocidade MENOR do que tiro.velocidade, com isso conseguiria    *
     * curvar e atingir na região de impactação com a inclinação medida.*
     * Só pode corrigir a altura depois que finalizar correção de theta.*
     ********************************************************************/
bool correcaoAltura(struct prjt *projetil, struct disparo *tiro, const double impacto_altura) {
    double delta_y = fabs (projetil->y - impacto_altura); //parâmetro para comparação entre a altura após atingir e a altura calculada após as iterações ao sair da edificação.
    if ( delta_y > 0.01){
        if (projetil->y > impacto_altura) {
            tiro->altura -= delta_y/2 + 0.001;   //delta_y >= 0 pelo "fabs" aplicado anteriormente.
            return true;
        } else {
            tiro->altura += delta_y/2 + 0.001;
            return true;
        }
    }
    return false;
}

void inicializarProjetilEW(double *t, struct prjt *projetil, const struct disparo *tiro, struct vento *w, const double yInicial) {
    *t=0;
    projetil->x = 0.0;
    projetil->y = yInicial; // Permite inicialização com um valor diferente de zero
    projetil->z = 0.0;
    projetil->vx = tiro->velocidade * cos(tiro->theta);
    projetil->vy = tiro->velocidade * sin(tiro->theta);
    projetil->vz = 0.0; // Velocidade inicial em desvio perpendicular ao downrange

    // Recalculando o vento pela variação do azimute inicial
    w->x = w->norte * cos(tiro->azimute) + w->leste * sin(tiro->azimute);
    w->y = 0.0;
    w->z = -w->norte * sin(tiro->azimute) + w->leste * cos(tiro->azimute);

    projetil->taxa_de_subida = atan2 (projetil->vy,projetil->vx);
    projetil->rumo = tiro->azimute + atan2 (projetil->vz,projetil->vx);

}


void inicializarTiro(const struct prjt *projetil, const struct impactacao *impacto, struct disparo *tiro, const struct vento *w, const struct edificacao *edificio, enum origem_disparo origem) {
    if (origem == Nivel_do_Mar){
        tiro->altura = 0;
        tiro->latitude = (180/M_PI)*impacto->latitude - distLatGraus (projetil, tiro);  //latitude ao Nível do Mar
        tiro->longitude = impacto->longitude - distLongGraus (projetil, impacto, tiro);  //longitude ao Nível do Mar
        tiro->origem = Nivel_do_Mar;
        // tiro->theta e tiro->azimute para o caso nível do mar é estimado na primeira iteração da função movimentoProjetil
    } else {
        tiro->altura = projetil->y;
        tiro->latitude = edificio->latitude;
        tiro->longitude = edificio->longitude;
        tiro->theta = projetil->taxa_de_subida;
        tiro->azimute = projetil->rumo;
        tiro->origem = Edificacao;
    }
}




//linha 233
double movimentoProjetil(int *n, struct prjt *projetil, struct impactacao *impacto, struct disparo *tiro, struct vento *w, struct edificacao *edificio, bool calcular_Edf){
 
    FILE *arquivo;
    arquivo = fopen("data","w");

    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        exit(1);
    }
    
    bool corrigirTheta = true;
    bool corrigirAzimute = true;
    bool corrigirAltura = true;
    bool ultimarodada = false;

    if(calcular_Edf == false) corrigirAltura = false; // só corrige altura se houver prédio

    double downrangeMax = projetil->x; // pra subtrair do Downrange do NMM até o prédio e ficar somente a distância entre o prédio e a impactacao. Valor só será diferente de 0 na segunda vez, ou seja, quando tiver edificio
    double distanciaPredio_Impacataco = 0;

    double t;
    if (calcular_Edf){
        lerDadosEdificacao("edificacao.txt", edificio);
        distanciaPredio_Impacataco = haversine (180*impacto->latitude/M_PI, impacto->longitude, edificio->latitude, edificio->longitude);
    }

    tiro->altura = 0; // incialmente, considera ao nível do mar, se hovuer edf: muda no decorrer

    while (corrigirTheta || corrigirAzimute || corrigirAltura || ultimarodada) {

        //Condições iniciais para variáveis com edificação no caminho:
        inicializarProjetilEW(&t, projetil, tiro, w, tiro->altura);
        if(ultimarodada) fprintf(arquivo,"A posição e altura do projétil dadas a partir do momento do disparo na seguinte ordem:"
                                            "\nTempo\tLatitude\t\t\tLongitude\t\tAltura:\n");
        /********************************************************************************************************************
         * Início do laço para cálculo RK4 para um disparo ocorrido de uma edificação mais próxima da região de impactação  *
         *                                                                                                                  *
         ********************************************************************************************************************/

        //while ( projetil.x < downrangeMax){ <- downrangeMax != "caminho total percorrido no ar"
        while ( (impacto->phi<0) ?   (fabs(projetil->y-impacto->altura)>0.1 || projetil->taxa_de_subida > 0) :   (projetil->y < impacto->altura) ){

            t += H;
            atualizarProjetil(projetil, w, calcularg(impacto->latitude), tiro->azimute);
            
            if (ultimarodada) {
                fprintf(arquivo, "%.3lf\t%.12lf, %.12lf\t %lf m\n", t, tiro->latitude + distLatGraus (projetil, tiro), tiro->longitude + distLongGraus (projetil, impacto, tiro), projetil->y);
                fflush(arquivo);
                }
            //Na equação das velocidades, uma das variáveis é o Azimute atual em relação ao Norte, por isso o termo recebe somado a "inclinação lateral", pois assim será o azimute naquela posição do projétil.
            //Se com arrasto, o projétil caiu antes de atingir a "altura" é porque para o dado theta, ele não subirá muito. Então deve-se incrementar "theta".
            if (projetil->y<0 ) {
                tiro->theta = tiro->theta + PARADA;
                if (tiro->theta>0.785) { //Se theta estiver a 45 graus e ainda sim o projétil não atingir a altura, não há solução para o problema com os parâmetros fornecidos. O programa encerra-se.
                    printf ("\nCom os parâmetros fornecidos, o projétil não atingiria a altura de impactação.");
                    printf("O programa não continuará a ser executado\n");
                    exit(1);
                }
                break; // necessário sair do while caso tiro->theta = tiro->theta + PARADA, para recálculo de inicialização de variáveis
            }

            if (calcular_Edf && projetil->x > (downrangeMax-distanciaPredio_Impacataco) ) {
            //essa condicao aumenta o erro se o usuário colocar uma edificação fora da trajetória
            //Essa condicao testa para ver se o projétil já passou da edificação. Nesse caso, testa (abaixo) para saber se o projétil passou por cima da edificação.

                if(projetil->y < edificio->altura){

                    t=0.0;
                    inicializarTiro(projetil, impacto, tiro, w, edificio, Edificacao);     
                    downrangeMax = downrangeMax - projetil->x; // subtraindo a posição atual do projétil, temos a distância da edificação até a impactacao.
                    inicializarProjetilEW(&t, projetil, tiro, w, tiro->altura);
                    calcular_Edf = false;                       // == 0 Interrompe os testes para esses cálculos (condição do "if" acima).


                } else {
                    printf ("O projétil passou por cima da edificação.");
                    inicializarTiro(projetil, impacto, tiro, w, edificio, Nivel_do_Mar);
                    inicializarProjetilEW(&t, projetil, tiro, w, tiro->altura);
                    calcular_Edf = false;
                }
            }
        }

        (*n)++;
        projetil->rumo = tiro->azimute + atan2 (projetil->vz,projetil->vx);

        // Etapa de correções de ângulo e altura inicial.
        if (corrigirTheta) {
            corrigirTheta = correcaoTheta(projetil, tiro, impacto->phi);
        }

        if (corrigirAltura && !corrigirTheta) {
            corrigirAltura = correcaoAltura(projetil, tiro, impacto->altura);
        }

        /*Checagem para ver se o projétil parte acima da edificação*/
        if (tiro->altura > edificio->altura + 1.5) { /*+1.5 porque uma pessoa, em tese, pode subir no topo e disparar. Cuidado porque o tiro pode vir de antes do prédio. Analisar com calma.*/
            printf("\nATENÇÃO!\n");
            printf("\nATENÇÃO!\n");
            printf("\nPara que o projétil atinja uma altura de %.2lf m com inclinação de %.0lfº, disparado a %.0lf m de distância da impactação (Coordenadas: %lf ,%lf), precisaria ter sido disparado a uma altura de %.2lf m do solo (considerando que sua velocidade inicial é %.0lf m/s), ou seja, mais alto que a suposta edificação de origem, logo, não pode ter partido desta edificação.",impacto->altura,180*impacto->phi/M_PI,distanciaPredio_Impacataco,edificio->latitude,edificio->longitude,tiro->altura,tiro->velocidade);
            printf("\n\nReveja as condições inciais e reinicie o programa.\n\n");
            printf("\nATENÇÃO!\n");
            printf("O programa não continuará a ser executado\n");
            exit(1);
        }

        /********************************************************************************************************
         * Etapa de correção do AZ0 (Azimute Inicial = tiro.azimute)                                            *
         *                                                                                                      *
         * Considerando que impacto.azimute (γ) é a inclinação lateral medida (em relação ao Norte em sentido   *
         * horário), tal qual o φ para a angulação vertical, para descobrir o azimute inicial do disparo, é     *
         * necessário que o rumo do projétil (projetil.rumo) - γ seja menor que um delta arbitrariamente        *
         * escolhido. Ou seja, como γ (tiro.azimute) é fixo (medido) e o rumo do projétil varia com a simulação,*
         * vamos corrigir o azimute (incrementando ou subtraindo) para obter um resultado menor que um delta    *
         * escolhido.                                                                                           *
         *                                                                                                      *
         ********************************************************************************************************/
        if (corrigirAzimute && !corrigirTheta) {
            corrigirAzimute = correcaoAzimute(projetil, tiro, impacto->azimute);
        }

        if (!corrigirTheta && !corrigirAzimute && !corrigirAltura){
            if (ultimarodada) break;
            ultimarodada = true;
        }
    }
    fclose(arquivo);
    return t;
}

double calcularg(double latitude){
    return 9.780327*(1+0.0053024*sin(latitude)*sin(latitude) - 0.0000058*sin(2*latitude)*sin(2*latitude));
}