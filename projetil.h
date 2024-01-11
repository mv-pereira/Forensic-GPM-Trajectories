#ifndef PROJETIL_H
#define PROJETIL_H

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PARADA 0.00174533    // Critério de parada para ajuste de angulação. 0.00174533 rad = 0.1º.
#define PARADA_ALTURA 0.01

#include <string.h>

// Estrutura do Projétil
enum sentido_rotacao {Dextrogiro, Levogiro};

struct caracteristicas_do_projetil {
    enum sentido_rotacao rotacao;
    double massa;
    double diametro;
    double coef_arrasto;
};

struct prjt {
    double x, y, z;
    double vx, vy, vz;
    double taxa_de_subida, rumo; // Inclinação e inclinação lateral instantânea
    double latitude, longitude;
    struct caracteristicas_do_projetil propriedades;
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
double distLatGraus(const struct prjt *projetil, const struct disparo *tiro);
double distLongGraus(const struct prjt *projetil, const struct impactacao *impacto, const struct disparo *tiro);
double spindrift(double tempo);
double pos_x(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_y(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double pos_z(double kappa, struct prjt *projetil, double inclinacao_RK_anterior, struct vento *w, double g);
double w_vx(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double w_vy(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double w_vz(double k, struct prjt *projetil, double correcao, struct vento *w, double g);
double runge_kutta(double (*funcao)(double, struct prjt *, double, struct vento *, double), struct prjt *projetil, struct vento *w, double passo, double kappa, double g);
double ajustar_theta(double phi_final, double phi_medido, double theta);
double ajuste_AZ(double azimute_disparo, double azimute_final, double azimete_medido);

#endif // PROJETIL_H
