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
    return 1.0/(cos (alpha));
}
double arcsec(double x){
    return acos(1.0/x); //arcsec t = arccos(1/t).
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
    double raio = raioTerraLatitude ((lat1+lat2)/2); //Raio da terra.
    return 2*raio*asin(sqrt( pow( sin((lat2-lat1)/2),2) + cos(lat1)*cos(lat2)*pow(sin((long2-long1)/2),2)));
}

// Este cálculo assume um modelo esférico da Terra, que é uma aproximação. 
// Para maior precisão, especialmente em grandes distâncias, um modelo elipsoidal seria mais adequado.
void navegEsferica(struct prjt *projetil, double distancia) {
    double angular_distance = distancia / raioTerraLatitude(projetil->latitude);
    double novaLat = asin(sin(projetil->latitude) * cos(angular_distance) + cos(projetil->latitude) * sin(angular_distance) * cos(projetil->rumo));
    double novaLong = projetil->longitude + atan2(sin(projetil->rumo) * sin(angular_distance) * cos(projetil->latitude), cos(angular_distance) - sin(projetil->latitude) * sin(novaLat));

    projetil->latitude = novaLat;
    projetil->longitude = novaLong;
}

/********************************************************************
 * Cáluclo da distância percorrida pelo projétil convertida         *
 * para graus latitude/longitude tendo em vista sua localização     *
 * relativa às edificações e impactações.                           *
 *                                                                  *
 * "distLat" é a "distância" em radianos entre a latitude do        *
 * disparo e a latitude do ponto solicitado.                        *
 * "distLong" é a "distância" em radianos entre a longitude do      *
 * disparo e a longitude do ponto solicitado.                       *
 *                                                                  *
 * Constante tem unidades de rad/m.                                 *
 ********************************************************************/

double distLat (const struct prjt *projetil, const struct disparo *tiro){
    double raio = raioTerraLatitude(projetil->latitude);
    double c = 360.0/(2*M_PI*raio);
    return c*(projetil->x*cos(tiro->azimute) - projetil->z*sin(tiro->azimute));
}

double distLong(const struct prjt *projetil, const struct impactacao *impacto, const struct disparo *tiro){
    double raio = raioTerraLatitude(projetil->latitude);
    double c = 360.0/(2*M_PI*raio);
    return c*cos(impacto->latitude)*(projetil->x*sin(tiro->azimute) + projetil->z*cos(tiro->azimute) );
}

/****************************************************************************
 *         SG - Miller Stability Formula: Velocity Correction               *
 * Fórmula original (v/2800)^(1/3), v em fps                                *
 * Applied Ballistics for Long Range Shooting, 3ed, pg 429, Appendix B      *
 ****************************************************************************/
double v_correction_msf(double v){
    return powl(v/853.44,1.0/3.0);
}
/****************************************************************************
 *                        SG - Miller Stability Formula                     *
 * Valores estimados para projéteis subsônicos a partir de                  *
 * Applied Ballistics for Long Range Shooting, 3ed, pg 428, Appendix B      *
 * twist -> Calibers/turn                                                   *
 * length -> bullet length in calibers                                      *
 ****************************************************************************/
double miller_stability_formula (const double mass, const double twist, const double diameter, const double length){    
    double m = 15.4324*mass;    // in grain
    double d = diameter/25.4;   // in inches
    double t = twist/d;         // in calibers/turn
    double l = (length/25.4)/d; // length in calibers
    return 30.0*m/(pow(t,2)*pow(d,3)*l*(1+l*l));
}

/****************************************************************************
 *                        Spindrift aproximado                              *
 * Valores estimados para projéteis subsônicos a partir de                  *
 * Applied Ballistics for Long Range Shooting, 3ed, pg 423, Appendix B      *
 ****************************************************************************/
double spindrift (double tempo, double sg){
    double drift;
    drift = 1.25*(sg+1.2)*pow(tempo,1.83);  //Resultado em polegada.
    return 0.0254*drift;                      //Resultado em metro.
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
    return (  -k*(sqrtl ( powl((projetil->vx +correcao -w->x),2) + powl((projetil->vy +correcao -w->y),2) + powl((projetil->vz +correcao -w->z),2) ))*(projetil->vz +correcao -w->z) + 2*OMEGA*( (projetil->vx +correcao -w->x)*(sin (projetil->latitude)) -(projetil->vy +correcao -w->y)*(cos (projetil->latitude))*(cos (projetil->rumo)) ));
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

/********************************************************************************************************
 * funções auxiliares para cálculo de posição com contribuição exclusiva do efeito coriolis.            *
 * Note que a função retorna apenas a parte relativa ao efeito coriolis das funções auxiliares w        *
 ********************************************************************************************************/

double pos_x_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return 2*OMEGA*( -(projetil->vy +inclinacao_RK_anterior -w->y)*(cos (projetil->latitude))*(sin (projetil->rumo)) -(projetil->vz +inclinacao_RK_anterior -w->z)*(sin (projetil->latitude)));
}

double pos_y_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return 2*OMEGA*( (projetil->vx +inclinacao_RK_anterior -w->x)*(cos (projetil->latitude))*(sin (projetil->rumo)) +(projetil->vz +inclinacao_RK_anterior -w->z)*(cos (projetil->latitude))*(sin (projetil->rumo)));
}

double pos_z_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g){
    return 2*OMEGA*( (projetil->vx +inclinacao_RK_anterior -w->x)*(sin (projetil->latitude)) -(projetil->vy +inclinacao_RK_anterior -w->y)*(cos (projetil->latitude))*(cos (projetil->rumo)));
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


double calcularg(double latitude){
    //https://pt.wikipedia.org/wiki/Acelera%C3%A7%C3%A3o_da_gravidade
    return 9.780327*(1+0.0053024*sin(latitude)*sin(latitude) - 0.0000058*sin(2*latitude)*sin(2*latitude));
}


void atualizarProjetil(struct prjt *projetil, struct vento *w, double g, double azimuteTiro, double d_spin) {
    //Constante de arrasto/ro. Note que a massa já entra na constante. Densidade do ar será calculado nas funcoes de vx,vy,vz...
    //kappa - Variável da qual depende o arrasto. Há o cálculo inicial de tal kappa sem densidade e a cada altura (projetil.y) é calculada a densidade do ar.
    double kappaSdensidade = (projetil->propriedades.coef_arrasto/projetil->propriedades.massa)*(projetil->propriedades.diametro/2); //Ordem alterada de c*a/(2.0*m) para evitar Underflow.
    double kappa = kappaSdensidade * densidade_ar(projetil->y); // kappaSdensidade deve ser definido em algum lugar do seu código

    if (projetil->cor.calcular) {
        projetil->cor.x += runge_kutta(&pos_x_cor, projetil, w, H, kappa, g) * H;
        projetil->cor.y += runge_kutta(&pos_y_cor, projetil, w, H, kappa, g) * H;
        projetil->cor.z += runge_kutta(&pos_z_cor, projetil, w, H, kappa, g) * H;
    }

    projetil->x += runge_kutta(&pos_x, projetil, w, H, kappa, g) * H;
    projetil->y += runge_kutta(&pos_y, projetil, w, H, kappa, g) * H;
    projetil->z += runge_kutta(&pos_z, projetil, w, H, kappa, g) * H + d_spin; //Foi calculado inicialmente o spin total e dividido pelo passo, para integrar no total.

    projetil->vx += runge_kutta(&w_vx, projetil, w, H, kappa, g) * H;
    projetil->vy += runge_kutta(&w_vy, projetil, w, H, kappa, g) * H;
    projetil->vz += runge_kutta(&w_vz, projetil, w, H, kappa, g) * H;

    projetil->taxa_de_subida = atan2(projetil->vy, projetil->vx);
    projetil->rumo = azimuteTiro + atan2(projetil->vz, projetil->vx);
    
    double distancia_horz_percorrida = sqrtl(pow(projetil->x,2)+pow(projetil->z,2));
    navegEsferica(projetil, distancia_horz_percorrida);
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
    if ( delta_y > PARADA_ALTURA){
        if (projetil->y > impacto_altura) {
            tiro->altura -= delta_y + PARADA_ALTURA/2;   //delta_y >= 0 pelo "fabs" aplicado anteriormente.
                if (tiro->altura < 0) {
                    tiro->altura = 0;
                return false;
                }
        } else {
            tiro->altura += delta_y + PARADA_ALTURA/2;
        }
        return true;
    }
    return false;
}


void inicializarProjetilEW(double *t, struct prjt *projetil, const struct disparo *tiro, struct vento *w, const double yInicial) {
    *t=0;
    if (projetil->cor.calcular) {
        projetil->cor.x = 0;
        projetil->cor.y = 0;
        projetil->cor.z = 0;
    }
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
    projetil->latitude = tiro->latitude;
    projetil->longitude = tiro->longitude;

}


void inicializarTiro(const struct prjt *projetil, const struct impactacao *impacto, struct disparo *tiro, const struct vento *w, const struct edificacao *edificio, enum origem_disparo origem) {
    if (origem == Nivel_do_Mar){
        tiro->altura = 0;
        tiro->latitude = impacto->latitude - distLat (projetil, tiro);  //latitude ao Nível do Mar
        tiro->longitude = impacto->longitude - distLong (projetil, impacto, tiro);  //longitude ao Nível do Mar
        tiro->origem = Nivel_do_Mar;
        // tiro->theta e tiro->azimute para o caso nível do mar é estimado na primeira iteração/leitura do arquivo
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
double movimentoProjetil(int *n, struct prjt *projetil, struct impactacao *impacto, struct disparo *tiro, struct vento *w, struct edificacao *edificio, bool calcular_Edf, double *downrangeMax, const double distanciaPredio_Impacataco, double d_spin){
 
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

    double t;
    

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
        while ( (impacto->phi<0) ?   (fabs(projetil->y-impacto->altura)>PARADA_ALTURA || projetil->taxa_de_subida > 0) :   (projetil->y < impacto->altura) ){

            t += H;
            atualizarProjetil(projetil, w, calcularg(projetil->latitude), tiro->azimute, d_spin);
            
            if (ultimarodada) {
                fprintf(arquivo, "%.3lf\t%.12lf, %.12lf\t %lf m\n", t, tiro->latitude + distLat (projetil, tiro), tiro->longitude + distLong (projetil, impacto, tiro), projetil->y);
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

            if (calcular_Edf && projetil->x > (*downrangeMax-distanciaPredio_Impacataco) ) {
            //essa condicao aumenta o erro se o usuário colocar uma edificação fora da trajetória
            //Essa condicao testa para ver se o projétil já passou da edificação. Nesse caso, testa (abaixo) para saber se o projétil passou por cima da edificação.

                if(projetil->y < edificio->altura){

                    t=0.0;
                    inicializarTiro(projetil, impacto, tiro, w, edificio, Edificacao);     
                    *downrangeMax = *downrangeMax - projetil->x; // subtraindo a posição atual do projétil, temos a distância da edificação até a impactacao.
                    inicializarProjetilEW(&t, projetil, tiro, w, tiro->altura);
                    calcular_Edf = false;                       // == 0 Interrompe os testes para esses cálculos (condição do "if" acima).


                } else {
                    printf ("O projétil passou por cima desta edificação.\n\n");
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
            printf("\nPara que o projétil atinja uma altura de %.2lf m com inclinação de %.0lfº, disparado a %.0lf m de distância da impactação (Coordenadas: %lf ,%lf), precisaria ter sido disparado a uma altura de %.2lf m do solo (considerando que sua velocidade inicial é %.0lf m/s), ou seja, mais alto que a suposta edificação de origem, logo, não pode ter partido desta edificação.",impacto->altura,180*impacto->phi/M_PI,distanciaPredio_Impacataco, edificio->latitude, edificio->longitude,tiro->altura,tiro->velocidade);
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
    if(calcular_Edf == false) *downrangeMax = projetil->x; // O cálculo terminava na primeira iteração e downrangeMax não era avaliado, pois não havia prédio
    fclose(arquivo);
    return t;
}