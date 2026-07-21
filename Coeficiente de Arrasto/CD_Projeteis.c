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
  * Programa para análise de curvas de arrasto Gx        *
  * aplicadas a projéteis de dimensões conhecidas.       *
  *                                                      *
  * O programa utiliza velocidades fornecidas pelo       *
  * fabricante em 0 m, 50 m e 100 m para comparar        *
  * diferentes curvas de arrasto, ajustar um fator       *
  * multiplicativo de Cd e selecionar a curva com menor  *
  * erro RMS.                                           *
  *                                                      *
  * Também permite exportar tabela Cd x velocidade,      *
  * útil para uso em outros programas balísticos.        *
  *                                                      *
  * Mario Pereira                                       *
  * Perito Criminal - DEPP/ICPAS - Pernambuco           *
  ********************************************************/

 #include <stdio.h>

 #include <math.h>

 #include <locale.h>

 #include <stdlib.h>

 #include <float.h>

 #include "curvas_gx.h"


 #define OMEGA 0.000072921 //Taxa de rotação da terra em "rad/s".
 #define PARADA 0.0174533 // Critério de parada para ajuste de angulação. 0.0174533 rad = 1º.
 #define H 0.0001 //passo da iteração do Runge-Kutta.
 #define DEBUG 0 //DEBUG == 1 permite iniciar o programa com parâmetros pré setados.

 #define DISTANCIA_50 50.0
 #define DISTANCIA_100 100.0

 #define TEMP_PADRAO_C 15.0

 #define FATOR_INICIAL 1.0
 #define FATOR_PASSO 0.01
 #define FATOR_TOLERANCIA 0.0001

 #define FATOR_MIN 0.10
 #define FATOR_MAX 5.00

 //Estrutura do Projétil

 enum sentido_rotacao {
   Dextrogiro,
   Levogiro
 };

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
   double x, y, z; //Componentes x,y,z na direção de deslocamento principal (downrange) do projétil.
   double norte, leste; //Componentes Norte e Leste da direção do vento.
 };

 //Estrutura da Edificação
 struct edificacao {
   double latitude;
   double longitude;
   double altura;
 };

 //Estrutura do Disparo
 enum origem_disparo {
   Nivel_do_Mar,
   Edificacao
 };

 struct disparo {
   enum origem_disparo origem;
   double latitude;
   double longitude;
   double altura;
   double velocidade; //Estimado pela característica do projétil.
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

 struct resultado_curva {
   const char * nome_curva;
   double fator_ajuste;
   double v1_calculada;
   double v2_calculada;
   double erro_1;
   double erro_2;
   double erro_rms;
 };

 /****************************************************************************
  *Aproximacao exponencial para Densidade                                    *
  *https://en.wikipedia.org/wiki/Density_of_air#Exponential_approximation    *
  ****************************************************************************/

 double densidade_ar(double altura) {
   return 1.22501236323 * exp(-altura / 10400);
 }

 double velocidade_som(double temperatura_celsius) {
   return 331.3 + 0.606 * temperatura_celsius;
 }

 double velocidade_relativa(struct prjt * projetil, struct vento * w) {
   return sqrt(
     pow(projetil -> vx - w -> x, 2) +
     pow(projetil -> vy - w -> y, 2) +
     pow(projetil -> vz - w -> z, 2)
   );
 }

 double calcula_mach(struct prjt * projetil, struct vento * w, double velocidade_do_som) {
   return velocidade_relativa(projetil, w) / velocidade_do_som;
 }

 /****************************************
  * Funções trigonométricas auxiliares   *
  ****************************************/

 double sec(double alpha) {
   return 1 / (cos(alpha));
 }
 double arcsec(double x) {
   return acos(1 / x); //arcsec t = arccos(1/t).
 }

 /****************************************************************
  * funções auxiliares para cálculo de posição no RK.            *
  ****************************************************************/

 double pos_x(double kappa, struct prjt * projetil, double inclinacao_RK_anterior, struct vento * w, double g) {
   return projetil -> vx;
 }

 double pos_y(double kappa, struct prjt * projetil, double inclinacao_RK_anterior, struct vento * w, double g) {
   return projetil -> vy;
 }

 double pos_z(double kappa, struct prjt * projetil, double inclinacao_RK_anterior, struct vento * w, double g) {
   return projetil -> vz;
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
 double w_vx(double k, struct prjt * projetil, double correcao, struct vento * w, double g) {
   return (-k * (sqrtl(powl((projetil -> vx + correcao - w -> x), 2) + powl((projetil -> vy + correcao - w -> y), 2) + powl((projetil -> vz + correcao - w -> z), 2))) * (projetil -> vx + correcao - w -> x) + 2 * OMEGA * (-(projetil -> vy + correcao - w -> y) * (cos(projetil -> latitude)) * (sin(projetil -> rumo)) - (projetil -> vz + correcao - w -> z) * (sin(projetil -> latitude))));
 }

 /*****************************************************************
  * funcao Auxiliar w_vy para Vel y (Altura): Esta função apenas  *
  *                                           calculra o 'k' para *
  *                                           projetil_1.vy.      *
  *                                                               *
  *****************************************************************/
 double w_vy(double k, struct prjt * projetil, double correcao, struct vento * w, double g) {
   return (-k * (sqrtl(powl((projetil -> vx + correcao - w -> x), 2) + powl((projetil -> vy + correcao - w -> y), 2) + powl((projetil -> vz + correcao - w -> z), 2))) * (projetil -> vy + correcao - w -> y) - g + 2 * OMEGA * ((projetil -> vx + correcao - w -> x) * (cos(projetil -> latitude)) * (sin(projetil -> rumo)) + (projetil -> vz + correcao - w -> z) * (cos(projetil -> latitude)) * (sin(projetil -> rumo))));
 }

 /******************************************************************
  * funcao Auxiliar para Vel z (Deriva/drift): Esta função apenas  *
  *                                            calculra o 'k' para *
  *                                            projetil_1.vz1.     *
  *                                                                *
  ******************************************************************/
 double w_vz(double k, struct prjt * projetil, double correcao, struct vento * w, double g) {
   return (-k * (sqrtl(powl((projetil -> vx + correcao - w -> x), 2) + powl((projetil -> vy + correcao - w -> y), 2) + powl((projetil -> vz + correcao - w -> z), 2))) * (projetil -> vz + correcao - w -> z) + 2 * OMEGA * ((projetil -> vx + correcao - w -> x) * (sin(projetil -> latitude)) - (projetil -> vy + correcao - w -> y) * (cos(projetil -> latitude)) * (cos(projetil -> rumo))));
 }

 /************************************************
  * Função para cálculo de Runge-Kutta 4a Ordem  *
  ************************************************/
 double runge_kutta(double( * funcao)(double, struct prjt( * ), double, struct vento( * ), double), struct prjt * projetil, struct vento * w, double passo, double kappa, double g) {
   double k, k1, k2, k3, k4;

   k1 = funcao(kappa, projetil, 0, w, g);
   k2 = funcao(kappa, projetil, k1 * (passo / 2), w, g);
   k3 = funcao(kappa, projetil, k2 * (passo / 2), w, g);
   k4 = funcao(kappa, projetil, k3 * (passo / 2), w, g);
   k = (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
   return k;
 }

 double simula_velocidade_final(
   const curva_arrasto_t * curva,
     double fator_ajuste,
     double distancia_alvo,
     struct disparo tiro,
     struct prjt projetil_base,
     struct vento w,
     double g,
     double velocidade_do_som
 ) {
   double t = 0.0;
   double kappa;
   double kappaSdensidade;
   double mach;
   double cd_base;
   double cd_ajustado;

   struct prjt projetil = projetil_base;
   projetil.latitude = tiro.latitude;

   projetil.x = 0.0;
   projetil.y = 1.0;
   projetil.z = 0.0;

   projetil.vx = tiro.velocidade;
   projetil.vy = 0.0;
   projetil.vz = 0.0;

   projetil.taxa_de_subida = atan2(projetil.vy, projetil.vx);
   projetil.rumo = tiro.azimute + atan2(projetil.vz, projetil.vx);

   while (projetil.x < distancia_alvo) {
     t += H;

     projetil.x += runge_kutta( & pos_x, & projetil, & w, H, 0, 0) * H;
     projetil.y += runge_kutta( & pos_y, & projetil, & w, H, 0, 0) * H;
     projetil.z += runge_kutta( & pos_z, & projetil, & w, H, 0, 0) * H;

     projetil.rumo = tiro.azimute + atan2(projetil.vz, projetil.vx);

     mach = calcula_mach( & projetil, & w, velocidade_do_som);

     cd_base = interpola_cd(curva, mach);
     cd_ajustado = cd_base * fator_ajuste;

     kappaSdensidade =
       (cd_ajustado / projetil.propriedades.massa) *
       (projetil.propriedades.diametro / 2.0);

     kappa = kappaSdensidade * densidade_ar(projetil.y);

     projetil.vx += runge_kutta( & w_vx, & projetil, & w, H, kappa, g) * H;
     projetil.vy += runge_kutta( & w_vy, & projetil, & w, H, kappa, g) * H;
     projetil.vz += runge_kutta( & w_vz, & projetil, & w, H, kappa, g) * H;

     projetil.taxa_de_subida = atan2(projetil.vy, projetil.vx);
   }

   return velocidade_relativa( & projetil, & w);
 }

 double calcula_erro_rms(
   const curva_arrasto_t * curva,
     double fator_ajuste,
     struct disparo tiro,
     struct prjt projetil,
     struct vento w,
     double g,
     double velocidade_do_som,
     double distancia_1,
     double v1_fabricante,
     double distancia_2,
     double v2_fabricante,
     double * v1_calculada,
     double * v2_calculada
 ) {
   double erro_1;
   double erro_2;

   * v1_calculada = simula_velocidade_final(
     curva,
     fator_ajuste,
     distancia_1,
     tiro,
     projetil,
     w,
     g,
     velocidade_do_som
   );

   * v2_calculada = simula_velocidade_final(
     curva,
     fator_ajuste,
     distancia_2,
     tiro,
     projetil,
     w,
     g,
     velocidade_do_som
   );

   erro_1 = * v1_calculada - v1_fabricante;
   erro_2 = * v2_calculada - v2_fabricante;

   return sqrt((erro_1 * erro_1 + erro_2 * erro_2) / 2.0);
 }

 struct resultado_curva ajusta_fator_curva(
   const curva_arrasto_t * curva,
     struct disparo tiro,
     struct prjt projetil,
     struct vento w,
     double g,
     double velocidade_do_som,
     double distancia_1,
     double v1_fabricante,
     double distancia_2,
     double v2_fabricante
 ) {
   struct resultado_curva resultado;

   double fator;
   double melhor_fator = FATOR_INICIAL;

   double erro;
   double menor_erro = DBL_MAX;

   double v1_calculada = 0.0;
   double v2_calculada = 0.0;

   double melhor_v1 = 0.0;
   double melhor_v2 = 0.0;

   for (fator = FATOR_MIN; fator <= FATOR_MAX; fator += FATOR_PASSO) {
     erro = calcula_erro_rms(
       curva,
       fator,
       tiro,
       projetil,
       w,
       g,
       velocidade_do_som,
       distancia_1,
       v1_fabricante,
       distancia_2,
       v2_fabricante, &
       v1_calculada, &
       v2_calculada
     );

     if (erro < menor_erro) {
       menor_erro = erro;
       melhor_fator = fator;
       melhor_v1 = v1_calculada;
       melhor_v2 = v2_calculada;
     }
   }

   resultado.nome_curva = curva -> nome;
   resultado.fator_ajuste = melhor_fator;

   resultado.v1_calculada = melhor_v1;
   resultado.v2_calculada = melhor_v2;

   resultado.erro_1 = melhor_v1 - v1_fabricante;
   resultado.erro_2 = melhor_v2 - v2_fabricante;

   resultado.erro_rms = menor_erro;

   return resultado;
 }

 void exporta_cd_velocidade_csv(
   const char * nome_arquivo,
     const curva_arrasto_t * curva,
       double fator_ajuste,
       double velocidade_inicial,
       double velocidade_final,
       double velocidade_do_som
 ) {
   FILE * arquivo;
   double v;
   double mach;
   double cd_base;
   double cd_ajustado;

   arquivo = fopen(nome_arquivo, "w");

   if (arquivo == NULL) {
     printf("\nErro: não foi possível criar o arquivo CSV: %s\n", nome_arquivo);
     return;
   }

   fprintf(
     arquivo,
     "velocidade_m_s,mach,cd_base,cd_ajustado,curva,fator_ajuste\n"
   );

   for (v = velocidade_inicial; v >= velocidade_final; v -= 1.0) {
     mach = v / velocidade_do_som;
     cd_base = interpola_cd(curva, mach);
     cd_ajustado = cd_base * fator_ajuste;

     fprintf(
       arquivo,
       "%.3f,%.6f,%.6f,%.6f,%s,%.6f\n",
       v,
       mach,
       cd_base,
       cd_ajustado,
       curva -> nome,
       fator_ajuste
     );
   }

   fclose(arquivo);

   printf("\nArquivo CSV exportado com sucesso: %s\n", nome_arquivo);
 }

 /****************************************************
  * Função Principal: Requer alguns iniciadores para *
  *                   calcular Cd, o Coeficiente     *
  *                   de Arrasto.                    *
  *                                                  *
  ****************************************************/

 int main() {
   setlocale(LC_ALL, "Portuguese"); //Utilizando caracteres e acentuação da língua portuguesa.

   printf("\nPrograma para cálculo do Coeficiente de Arrasto Cd (adimensional) para projéteis de dimensões conhecidas utilizando parâmetros fornecidos pelos fabricantes de projétis.\n");

   double distancia_1, distancia_2;
   double v1_fabricante, v2_fabricante;
   double temperatura_celsius;
   double velocidade_do_som;
   double g;
   size_t i;

   struct prjt projetil; //Estrutura do projétil.
   struct vento w; //Definição da struct do vento.
   struct impactacao impacto; //Estrutura da Impactação.
   struct disparo tiro; //Estrutura do Tiro.

   tiro.latitude = -23.698308; //Latitude da Fábrica da CBC em São Paulo
   tiro.azimute = 0.0;

   #if DEBUG //Valores de teste
   tiro.velocidade = 367.0;
   distancia_1 = 50.0;
   v1_fabricante = 331.0;
   distancia_2 = 100.0;
   v2_fabricante = 306.0;
   temperatura_celsius = TEMP_PADRAO_C;
   velocidade_do_som = velocidade_som(temperatura_celsius);

   projetil.propriedades.massa = 11.66 / 1000;
   projetil.propriedades.diametro = M_PI * powl((10 / 1000.0), 2) / 4;

   #else
   /********************************
    * Características do projetil  *
    ********************************/

   printf("\nDigite a velocidade inicial V0, na distância de 0 m (em m/s):\n");
   scanf("%lf", & tiro.velocidade);

   printf("\nDigite a primeira distância informada pelo fabricante (em m):\n");
   scanf("%lf", & distancia_1);

   printf("\nDigite a velocidade do projétil nessa primeira distância (em m/s):\n");
   scanf("%lf", & v1_fabricante);

   printf("\nDigite a segunda distância informada pelo fabricante (em m):\n");
   scanf("%lf", & distancia_2);

   printf("\nDigite a velocidade do projétil nessa segunda distância (em m/s):\n");
   scanf("%lf", & v2_fabricante);

   printf("\nDigite a temperatura ambiente (em graus Celsius). Use 15 para atmosfera padrão:\n");
   scanf("%lf", & temperatura_celsius);

   velocidade_do_som = velocidade_som(temperatura_celsius);

   printf("\nDigite a massa do projétil (em gramas):\n");
   scanf("%lf", & projetil.propriedades.massa);
   projetil.propriedades.massa = projetil.propriedades.massa / 1000.0;

   printf("\nDigite o diâmetro do projétil (em milímetros):\n");
   scanf("%lf", & projetil.propriedades.diametro);

   projetil.propriedades.diametro =
     M_PI * powl((projetil.propriedades.diametro / 1000.0), 2) / 4;

   #endif

   //Valor médio estimado para início dos cálculos.
   projetil.propriedades.coef_arrasto = 0.2;

   //Aceleração da gravidade na latitude. (em m/s^2)
   g = 9.780327 * (1 + 0.0053024 * sin(tiro.latitude) * sin(tiro.latitude) - 0.0000058 * sin(2 * tiro.latitude) * sin(2 * tiro.latitude));

   printf("\nCurvas de arrasto carregadas:\n");

   for (i = 0; i < n_curvas_arrasto; i++) {
     printf(
       "%zu - %s com %zu pontos\n",
       i + 1,
       curvas_arrasto[i].nome,
       curvas_arrasto[i].n_pontos
     );
   }

   w.velocidade = 0.0;
   w.direcao = 0.0;
   w.x = 0.0;
   w.y = 0.0;
   w.z = 0.0;
   w.norte = 0.0;
   w.leste = 0.0;

   struct resultado_curva resultados[16];

   printf("\nAjuste automático do fator para cada curva:\n");

   printf(
     "\n%-5s | %-14s | %-12s | %-15s | %-13s | %-16s | %-12s\n",
     "Curva",
     "Fator",
     "V1 calc",
     "Erro 1",
     "V2 calc",
     "Erro 2",
     "Erro RMS"
   );

   printf(
     "-------------------------------------------------------------------------------------------------------------\n"
   );

   for (i = 0; i < n_curvas_arrasto; i++) {

     resultados[i] = ajusta_fator_curva( &
       curvas_arrasto[i],
       tiro,
       projetil,
       w,
       g,
       velocidade_do_som,
       distancia_1,
       v1_fabricante,
       distancia_2,
       v2_fabricante
     );

     printf(
       "%-5s | fator = %.4f | V(%.0f m) = %.2f | erro1 = %+6.2f | V(%.0f m) = %.2f | erro2 = %+6.2f | RMS = %.4f m/s\n",
       resultados[i].nome_curva,
       resultados[i].fator_ajuste,
       distancia_1,
       resultados[i].v1_calculada,
       resultados[i].erro_1,
       distancia_2,
       resultados[i].v2_calculada,
       resultados[i].erro_2,
       resultados[i].erro_rms
     );
   }

   size_t indice_melhor = 0;
   double menor_erro_global = DBL_MAX;

   for (i = 0; i < n_curvas_arrasto; i++) {
     if (resultados[i].erro_rms < menor_erro_global) {
       menor_erro_global = resultados[i].erro_rms;
       indice_melhor = i;
     }
   }

   printf("\nMelhor curva encontrada:\n");
   printf("---------------------------------------------\n");

   printf("Curva: %s\n", resultados[indice_melhor].nome_curva);
   printf("Fator de ajuste: %.4f\n", resultados[indice_melhor].fator_ajuste);

   printf("Distância 1: %.2f m\n", distancia_1);
   printf("Velocidade fabricante na distância 1: %.2f m/s\n", v1_fabricante);
   printf("Velocidade calculada na distância 1:  %.2f m/s\n", resultados[indice_melhor].v1_calculada);
   printf("Erro na distância 1: %+.2f m/s\n", resultados[indice_melhor].erro_1);

   printf("Distância 2: %.2f m\n", distancia_2);
   printf("Velocidade fabricante na distância 2: %.2f m/s\n", v2_fabricante);
   printf("Velocidade calculada na distância 2:  %.2f m/s\n", resultados[indice_melhor].v2_calculada);
   printf("Erro na distância 2: %+.2f m/s\n", resultados[indice_melhor].erro_2);

   printf("Erro RMS: %.4f m/s\n", resultados[indice_melhor].erro_rms);

   printf(
     "\nObservação: o erro RMS foi calculado apenas com os dois pontos de referência informados pelo usuário.\n"
   );

   printf(
     "Uso recomendado preferencialmente entre %.2f m/s e %.2f m/s.\n",
     tiro.velocidade,
     fmin(v1_fabricante, v2_fabricante)
   );

   double velocidade_minima_tabela;

   velocidade_minima_tabela = fmin(v1_fabricante, v2_fabricante);

   exporta_cd_velocidade_csv(
     "cd_velocidade_melhor_curva.csv", &
     curvas_arrasto[indice_melhor],
     resultados[indice_melhor].fator_ajuste,
     tiro.velocidade,
     velocidade_minima_tabela,
     velocidade_do_som
   );

   return 0;
 }