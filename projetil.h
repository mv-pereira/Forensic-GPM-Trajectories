#ifndef PROJETIL_H
#define PROJETIL_H

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h> //Para função exit() na condicional da abertura do arquivo;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PARADA 0.00174533    // Critério de parada para ajuste de angulação. 0.00174533 rad = 0.1º.
#define PARADA_ALTURA 0.01
#define H 0.0001            //passo da iteração do Runge-Kutta.

// Estrutura do Projétil
enum sentido_rotacao {Dextrogiro = 1, Levogiro = -1};

struct caracteristicas_do_projetil {
    enum sentido_rotacao rotacao;
    double massa;
    double diametro;
    double coef_arrasto;
    double twist;
    double length;
};

struct coriolis {
    bool calcular;
    double x;
    double y;
    double z;
};

struct prjt {
    double x, y, z;
    double vx, vy, vz;
    double taxa_de_subida, rumo; // Inclinação e inclinação lateral instantânea
    double latitude, longitude;
    double sg;
    struct caracteristicas_do_projetil propriedades;
    struct coriolis cor;
};

// Estrutura do Vento
struct vento {
    double velocidade, direcao;
    double x, y, z;           // Componentes x, y, z na direção de deslocamento principal do projétil.
    double norte, leste;      // Componentes Norte e Leste da direção do vento.
};

// Estrutura da Edificação
struct edificacao {
    double latitude;
    double longitude;
    double altura;
};

// Estrutura do Disparo
enum origem_disparo {Nivel_do_Mar, Edificacao};

struct disparo {
    enum origem_disparo origem;
    double latitude;
    double longitude;
    double altura;
    double velocidade;      // Estimado pela característica do projétil.
    double azimute;
    double theta;
};

// Estrutura da Impactação
struct impactacao {
    double latitude;
    double longitude;
    double altura;
    double phi;
    double azimute;
};

// Protótipos de funções
double densidade_ar(double altura);
double sec(double alpha);
double arcsec(double x);
double haversine(double lat1, double long1, double lat2, double long2);
double raioTerraLatitude(double latitude);
double distLat(const struct prjt *projetil, const struct disparo *tiro);
double distLong(const struct prjt *projetil, const struct impactacao *impacto, const struct disparo *tiro);
double spindrift(double tempo, double sg);
double pos_x(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_y(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_z(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double w_vx(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double w_vy(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double w_vz(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double pos_x_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_y_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_z_cor (double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double runge_kutta(double (*funcao)(double, struct prjt *, double, struct vento *, double), struct prjt *projetil, struct vento *w, double passo, double kappa, double g);
double ajustar_theta(double phi_final, double phi_medido, double theta);
double ajuste_AZ(double azimute_disparo, double azimute_final, double azimete_medido);
double v_correction_msf(double v);
double miller_stability_formula (const double mass, const double twist, const double diameter, const double length);

// Calcula a aceleração da gravidade com base na latitude.
double calcularg(double latitude);
// Atualiza o estado do projétil com base nas condições atuais e no disparo.
void atualizarProjetil(struct prjt *projetil, struct vento *w, double g, double azimuteTiro, double d_spin);
bool correcaoTheta(struct prjt *projetil, struct disparo *tiro, const double impacto_phi);
bool correcaoAzimute(struct prjt *projetil, struct disparo *tiro, const double impacto_azimute);
bool correcaoAltura(struct prjt *projetil, struct disparo *tiro, const double impacto_altura);
// Inicializa o projétil e o vento com valores iniciais para simulação de disparo.
void inicializarProjetilEW(double *t, struct prjt *projetil, const struct disparo *tiro, struct vento *w, const double yInicial);
// Prepara o disparo com base nas informações do projétil, impactação, vento e edificação.
void inicializarTiro(const struct prjt *projetil, const struct impactacao *impacto, struct disparo *tiro, const struct vento *w, const struct edificacao *edificio, enum origem_disparo origem);
// Calcula o movimento do projétil e retorna informações relevantes sobre a trajetória.
double movimentoProjetil(int *n, struct prjt *projetil, struct impactacao *impacto, struct disparo *tiro, struct vento *w, struct edificacao *edificio, bool calcular_Edf, double *downrangeMax, const double distanciaPredio_Impacataco, double spin);

#endif // PROJETIL_H
