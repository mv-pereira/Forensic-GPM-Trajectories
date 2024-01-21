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

#define _GNU_SOURCE // Para usar qsort_r

#include <locale.h> //Utilizando caracteres e acentuação da língua portuguesa.
#include <time.h>   //Para velocidade de funções - DEBUG
#include "projetil.h"

struct listaEdificacoes {
    struct edificacao *edificios;
    double *distPredioImpact;
    int numEdificios;
};

void lerDadosDoArquivo(char *nomeArquivo, struct impactacao *impacto, struct disparo *tiro, struct prjt *projetil, struct vento *w, double *g);

void lerDadosEdificacao(char *nomeArquivo, struct edificacao *edificio);

// Calcula as distâncias entre as edificações e o ponto de impacto.
void calcularDistancias(struct listaEdificacoes *lista, struct impactacao *impacto);
void lerDadosEdificacao(char *nomeArquivo, struct edificacao *edificio);
struct listaEdificacoes lerDadosEdificacaoVar(char *nomeArquivo);

// Ordena as edificações com base nas distâncias do ponto de impacto.
void ordenarEdificacoes(struct listaEdificacoes *lista);

// Função de comparação usada por qsort_r para ordenar um array.
int compare(const void *a, const void *b, void *arg);

double pontoIntermediarioEDist(double lat1, double long1, double lat2, double long2, double lat3, double long3, double *latProx, double *longProx);


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

    char *nomeArquivo = "input";
    printf("\n------------------------------------------------------------------------------------------------------------\n");
    printf("Início da leitura dos dados a partir do arquivo: %s\n", nomeArquivo);
    printf("------------------------------------------------------------------------------------------------------------\n");
    lerDadosDoArquivo(nomeArquivo, &impacto, &tiro, &projetil, &w, &g);

/************************************
 * Início do cálculo SEM arrasto    *
 * para estimativa do θ inicial     *
 * a ser utilizado nas iterações    *
 ************************************/    
    printf("\n\n------------------------------------------------------------------------------------------------------------\n");
    printf("Cálculo do sistema sem arrasto.\n");
    printf("------------------------------------------------------------------------------------------------------------\n");


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

    printf("\n\n------------------------------------------------------------------------------------------------------------\n");
    printf("Leitura e ordenamento das edificações a partir da mais próxima da impactação.\n");
    printf("------------------------------------------------------------------------------------------------------------\n");
    struct listaEdificacoes edificacoes = lerDadosEdificacaoVar("edificacao");
    calcularDistancias(&edificacoes, &impacto);
    ordenarEdificacoes(&edificacoes);


    t = movimentoProjetil(&n, &projetil, &impacto, &tiro, &w, &edificacoes.edificios[0], calcular_Edf, &downrangeMax, edificacoes.distPredioImpact[0]);
    inicializarTiro(&projetil, &impacto, &tiro, &w, &edificio, Nivel_do_Mar);

    printf("Sem edificações e considerando o arrasto, os cálculos terminaram com os seguintes valores:"
            "\nTempo de deslocamento total do projétil: %.1f segundos."
            "\nDownrange Total = %.1lf m."
            "\nÂngulo θ (inicial) do disparo = %.1lfº."
            "\nÂngulo Az (inicial) = %.1lfº."
            "\nCoordenadas Georgráficas do disparo ao NMM (N/S, L/O): %.6lf, %.6lf.\n\n",
            t, projetil.x, (180/M_PI)*tiro.theta, (180/M_PI)*tiro.azimute, tiro.latitude, tiro.longitude);

    // Guardar as coordenadas do disparo partindo do solo para calcular a distância entre o edifício fornecido
    // e a trajetória calculada.
    double tiroLatSolo = tiro.latitude, tiroLongSolo = tiro.longitude;


/************************************
 * Fim das primeiras simulações     *
 * considerandoo projétil partindo  *
 * do solo.                         *
 ************************************/
    double latitudeMaisProxima, longitudeMaisProxima;
    if (edificacoes.numEdificios != 0) calcular_Edf = true;

    for (int i=0; i<edificacoes.numEdificios; i++){
        // Principal função

        printf("\n\n------------------------------------------------------------------------------------------------------------\n");
        printf("Testando a edificação #%d.\n",i+1);
        printf("Edifício %d: Coordenadas: (N/S, L/O): %.6lf, %.6lf, Altura %.2lf m, Distância até Impacto: %.2lf m.\n"
        "A distância lateral aproximada que o projétil partindo do solo (ou nível do mar) passa deste ponto vale: %.2f m.\n",
        i + 1,
        edificacoes.edificios[i].latitude*(180/M_PI),
        edificacoes.edificios[i].longitude*(180/M_PI),
        edificacoes.edificios[i].altura,
        edificacoes.distPredioImpact[i],
        pontoIntermediarioEDist(tiroLatSolo, tiroLongSolo, impacto.latitude, impacto.longitude, edificacoes.edificios[i].latitude, edificacoes.edificios[i].longitude, &latitudeMaisProxima, &longitudeMaisProxima));
        printf("------------------------------------------------------------------------------------------------------------\n");        
        t = movimentoProjetil(&n, &projetil, &impacto, &tiro, &w, &edificacoes.edificios[i], calcular_Edf, &downrangeMax, edificacoes.distPredioImpact[i]);
        
        if (tiro.origem == Edificacao) {
            printf ("O projétil partiu desta edificação\n\n");
            break; // Já tá calculando a partir do prédio mais próximo.
        }
    }
    // latitudeMaisP e longitudeMaisP foram criadas porque, ao se colocar manualmente a localização do edifício, era possível que o ponto marcado
    // não fosse o local por onde o projétil passaria, aumentando o erro a medida que a pessoa indica um ponto mais afastado do que o que o projétil passou
    // dessa forma, agora se está mostrando a rota do projétil, não exatamente onde a pessoa indicou. 
    if (tiro.origem == Nivel_do_Mar){
        latitudeMaisProxima = tiro.latitude;
        longitudeMaisProxima = tiro.longitude;
    }


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
           "\nAzimute inicial do disparo = %.2lfº.\n",n, projetil.x, projetil.y,(projetil.z<0 ? "esquerda" : "direita"), fabs(projetil.z), (180/M_PI)*tiro.theta, (180/M_PI)*projetil.taxa_de_subida, (180/M_PI)*tiro.azimute);
    printf("\nO efeito de Spin Drift foi responsável por deslocar o projétil aproximadamente, %.2f m para %s.\n", spindrift(t, projetil.sg), projetil.propriedades.rotacao == Dextrogiro ? "direta" : "esquerda");
    
    if (projetil.cor.calcular) {
        printf("\nO efeito Coriolis foi responsável por deslocar o projétil:\nEixo x: %f cm;\nEixo y: %f cm;\nEixo z: %f cm;\n",100*projetil.cor.x, 100*projetil.cor.y, 100*projetil.cor.z);
    } else{
        printf("\nO efeito Corilis não foi calculado separadamente.\n");
    }

    printf("\n\n------------------------------------------------------------------------------------------------------------\n");
    printf("O projétil partiu, aproximadamente, das coordenadas: %lf N/S, %lf L/O, a uma altura de %.2lf m, partindo %s.\n", (180/M_PI)*latitudeMaisProxima, (180/M_PI)*longitudeMaisProxima/*(180/M_PI)*tiro.latitude, (180/M_PI)*tiro.longitude*/, tiro.altura, (tiro.origem == Nivel_do_Mar ? "ao nível do mar" : "de uma edificação" ) ); //latitude_disparo e longitude_disparo foram calculados ao fim do segundo laço.
 
    velocidade_final = sqrt (pow(projetil.vx,2)+pow(projetil.vy,2)+pow(projetil.vz,2)); //Módulo nas três componentes.
    printf("\nO projétil levou cerca de %.1f segundos do momento do disparo à impactação.\nPossuía velocidade final de %.2f m/s e energia cinética de %.2f J.\n",t,velocidade_final,0.5*projetil.propriedades.massa*pow(velocidade_final,2));
    printf("O disparo foi efetuado com o cano apontado, aproximadamente, a %.0lf m acima do alvo, ou %.2lfº acima da impactação.\n", tan(tiro.theta)*projetil.x-impacto.altura, (180/M_PI)*(tiro.theta - atan2(impacto.altura, projetil.x)));
    printf("------------------------------------------------------------------------------------------------------------\n");    
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
    free(edificacoes.distPredioImpact);
    free(edificacoes.edificios);
    
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
    printf("Altura de Impacto: %.3lf m;\n", impacto->altura);

    fscanf(file, "%*[^:]: %lf", &impacto->phi);
    printf("Ângulo Phi: %.1lf graus;\n", impacto->phi);
    impacto->phi = impacto->phi*M_PI/180;

    fscanf(file, "%*[^:]: %lf", &impacto->azimute);
    printf("Ângulo Azimute: %.1lf graus;\n", impacto->azimute);
    impacto->azimute = impacto->azimute*M_PI/180;

    fscanf(file, "%*[^:]: %lf", &impacto->latitude);
    printf("Latitude do Impacto: %lf graus;\n", impacto->latitude);
    impacto->latitude = impacto->latitude*M_PI/180.0;

    *g = calcularg(impacto->latitude); //Açeleração da gravidade na latitude. (em m/s^2)

    fscanf(file, "%*[^:]: %lf", &impacto->longitude);
    printf("Longitude do Impacto: %lf graus;\n", impacto->longitude);
    impacto->longitude = impacto->longitude*M_PI/180;

    tiro->azimute = impacto->azimute;
    // Primeira estimativa da latitude e longitude do projétil é igualá-la ao da impactação
    // Depois de um ciclo, será atualizado, inicialmente, para a região do disparo
    projetil->latitude = impacto->latitude;
    projetil->longitude = impacto->longitude;
    tiro->latitude = impacto->latitude;
    tiro->longitude = impacto->longitude;

    fscanf(file, "%*[^:]: %lf", &tiro->velocidade);
    printf("Velocidade Inicial: %.2lf m/s;\n", tiro->velocidade);

    double massa_g;
    fscanf(file, "%*[^:]: %lf", &massa_g);
    printf("Massa do Projétil: %.3lf g;\n", massa_g);
    projetil->propriedades.massa = massa_g/1000.0;

    double diametro_mm;
    fscanf(file, "%*[^:]: %lf", &diametro_mm);
    printf("Diâmetro do Projétil: %.3lf mm;\n", diametro_mm);
    projetil->propriedades.diametro = M_PI*powl((diametro_mm/1000.0),2)/4;  // A divisão por quatro leva em conta o raio.

    fscanf(file, "%*[^:]: %lf", &projetil->propriedades.coef_arrasto);
    printf("Coeficiente de Arrasto: %lf;\n", projetil->propriedades.coef_arrasto);

    int rotacao;
    fscanf(file, "%*[^:]: %d", &rotacao);
    projetil->propriedades.rotacao = (rotacao == 1) ? Dextrogiro : Levogiro;
    printf("Rotação do Projétil: %d;\n", rotacao);

    fscanf(file, "%*[^:]: 1:%lf", &projetil->propriedades.twist);
    printf("Twist do Projétil: 1:%.1lf;\n", projetil->propriedades.twist);
    // projetil->propriedades.twist = (25.4*projetil->propriedades.twist)/diametro_mm;
    // printf("Twist do Projétil: %lf cal/turn;\n", projetil->propriedades.twist);

    fscanf(file, "%*[^:]: %lf", &projetil->propriedades.length);
    printf("Comprimento do Projétil: %.3lf mm;\n", projetil->propriedades.length);
    // projetil->propriedades.length = projetil->propriedades.length/diametro_mm; //length em calibers
    // printf("Comprimento do Projétil em unidades de Calibre: %lf;\n", projetil->propriedades.length);

    fscanf(file, "%*[^:]: %lf", &w->velocidade);
    w->velocidade = w->velocidade / 3.6; // Convertendo km/h para m/s
    printf("Velocidade do Vento: %.2lf m/s;\n", w->velocidade);

    fscanf(file, "%*[^:]: %lf", &w->direcao);
    printf("Direção do Vento: %.2lf graus;\n", w->direcao);
    w->direcao = (w->direcao+180.0)*M_PI/180.0;

    int coriolis;
    fscanf(file, "%*[^:]: %d", &coriolis);
    projetil->cor.calcular = (coriolis == 1) ? true : false;
    printf("O efeito Coriolis %s calculado separadamente.\n", (projetil->cor.calcular) ? "será" : "não será");

    projetil->sg = miller_stability_formula(massa_g, projetil->propriedades.twist, diametro_mm, projetil->propriedades.length)*v_correction_msf(tiro->velocidade);

    printf("Estabilidade Giroscópica estimada: %.2lf;\n",projetil->sg);
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

    edificio->latitude = edificio->latitude*(M_PI/180);
    edificio->longitude = edificio->longitude*(M_PI/180);

    fscanf(file, "%*[^:]: %lf", &edificio->altura);
    printf("Altura da Edificação: %lf metros\n", edificio->altura);

    fclose(file);
}

struct listaEdificacoes lerDadosEdificacaoVar(char *nomeArquivo) {
    FILE *file = fopen(nomeArquivo, "r");
    char buffer[256];
    char *token, *start;
    const char delim[2] = ",";
    struct listaEdificacoes lista;
    int i;

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

    // Lê o número de edificações
    //fgets(buffer, sizeof(buffer), file);
    sscanf(buffer, "Numero de edificacoes: %d", &lista.numEdificios);

    // Aloca memória para as edificações
    lista.edificios = (struct edificacao *)malloc(lista.numEdificios * sizeof(struct edificacao));
    if (lista.edificios == NULL) {
        printf("Erro na alocação de memória.\n");
        fclose(file);
        exit(1);
    }

    // Latitude
    fgets(buffer, sizeof(buffer), file);
    start = strchr(buffer, ':');
    if (start != NULL) {
        token = strtok(start + 1, delim);
        for (i = 0; i < lista.numEdificios; i++) {
            lista.edificios[i].latitude = (M_PI/180)*atof(token);
            token = strtok(NULL, delim);
        }
    }

    // Longitude
    fgets(buffer, sizeof(buffer), file);
    start = strchr(buffer, ':');
    if (start != NULL) {
        token = strtok(start + 1, delim);
        for (i = 0; i < lista.numEdificios; i++) {
            lista.edificios[i].longitude = (M_PI/180)*atof(token);
            token = strtok(NULL, delim);
        }
    }

    // Altura
    fgets(buffer, sizeof(buffer), file);
    start = strchr(buffer, ':');
    if (start != NULL) {
        token = strtok(start + 1, delim);
        for (i = 0; i < lista.numEdificios; i++) {
            lista.edificios[i].altura = atof(token);
            token = strtok(NULL, delim);
        }
    }

    fclose(file);
    return lista;
}

void calcularDistancias(struct listaEdificacoes *lista, struct impactacao *impacto) {
    lista->distPredioImpact = (double *)malloc(lista->numEdificios * sizeof(double));
    if (lista->distPredioImpact == NULL) {
        printf("Erro na alocação de memória.\n");
        exit(1);
    }

    for (int i = 0; i < lista->numEdificios; i++) {
        lista->distPredioImpact[i] = haversine(impacto->latitude, impacto->longitude, lista->edificios[i].latitude, lista->edificios[i].longitude);
    }
}

int compare(const void *a, const void *b, void *arg) {
    const double *distPredioImpact = (double *)arg;
    int indexA = *(const int *)a;
    int indexB = *(const int *)b;
    double distA = distPredioImpact[indexA];
    double distB = distPredioImpact[indexB];

    // Para ordenação decrescente
    if (distA > distB) return 1;
    if (distA < distB) return -1;
    return 0;
}

void ordenarEdificacoes(struct listaEdificacoes *lista) {
    int *indices = (int *)malloc(lista->numEdificios * sizeof(int));
    if (indices == NULL) {
        printf("Erro na alocação de memória.\n");
        exit(1);
    }

    // Inicializa os índices
    for (int i = 0; i < lista->numEdificios; i++) {
        indices[i] = i;
    }

    // Ordena os índices baseados nas distâncias
    qsort_r(indices, lista->numEdificios, sizeof(int), compare, lista->distPredioImpact);

    // Cria cópias temporárias para reordenar os arrays originais
    struct edificacao *edificiosOrdenados = (struct edificacao *)malloc(lista->numEdificios * sizeof(struct edificacao));
    double *distanciasOrdenadas = (double *)malloc(lista->numEdificios * sizeof(double));
    if (edificiosOrdenados == NULL || distanciasOrdenadas == NULL) {
        printf("Erro na alocação de memória.\n");
        exit(1);
    }

    // Reordena os arrays
    for (int i = 0; i < lista->numEdificios; i++) {
        edificiosOrdenados[i] = lista->edificios[indices[i]];
        distanciasOrdenadas[i] = lista->distPredioImpact[indices[i]];
    }

    // Copia os arrays ordenados de volta para os originais
    memcpy(lista->edificios, edificiosOrdenados, lista->numEdificios * sizeof(struct edificacao));
    memcpy(lista->distPredioImpact, distanciasOrdenadas, lista->numEdificios * sizeof(double));

    // Libera a memória temporária
    free(indices);
    free(edificiosOrdenados);
    free(distanciasOrdenadas);
}

// Para o cálculo da distância entre o ponto do edifício (escolhido pelo usuário) e a trajetória do projétil
// É necessário calcular a projeção do ponto C (edifício) na SemiReta do tiro (AB). O ponto projetado na
// SemiReta AB é D (latitude e longitude possível do disparo, já que estava na linha de tiro original). 
void pontoDprojecaoDeCemAB(double xa, double ya, double xb, double yb, double xc, double yc, double *xd, double *yd) {
    double abx, aby, acx, acy, t;

    // Vetores AB e AC
    abx = xb - xa;
    aby = yb - ya;
    acx = xc - xa;
    acy = yc - ya;

    // Calculando o parâmetro t para a projeção
    t = (acx * abx + acy * aby) / (abx * abx + aby * aby);

    // Garantindo que a projeção esteja na semirreta AB
    if (t < 0) t = 0;

    // Ponto D
    *xd = xa + t * abx;
    *yd = ya + t * aby;
}

double pontoIntermediarioEDist(double lat1, double long1, double lat2, double long2, double lat_edf, double long_edf, double *latProx, double *longProx) {
    pontoDprojecaoDeCemAB(long1, lat1, long2, lat2, long_edf, lat_edf, longProx, latProx);
    return haversine(*latProx, *longProx, lat_edf, long_edf);
}