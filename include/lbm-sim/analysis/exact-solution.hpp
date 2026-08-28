#ifndef __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP
#define __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP

#include "lbm-sim/analysis/types.hpp"

namespace lbm {
namespace analysis {

/**
 * Soluzione analitica per il Flusso di Couette (2D).
 * Parete inferiore fissa (y = 0), parete superiore mobile (y = H) con velocità
 * Umax. Profilo di velocità lineare: u_x(y) = Umax * (y / H), u_y = 0.
 */
class CouetteSolution2D : public Function<2> {
private:
  double H;    // Altezza del canale
  double Umax; // Velocità della parete mobile superiore

public:
  CouetteSolution2D(double channel_height, double max_velocity)
      : H(channel_height), Umax(max_velocity) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    // p.y rappresenta la coordinata verticale[cite: 3]
    double ux = Umax * (p.y / H);
    return utils::Vector<double, 2>{ux, 0.0};
  }
};

/**
 * Soluzione analitica per il Flusso di Poiseuille (2D).
 * Guidato da gradiente di pressione / forza di volume tra due pareti fisse (y =
 * 0 e y = H). Profilo di velocità parabolico: u_x(y) = 4 * Umax * (y / H) * (1
 * - y / H), u_y = 0.
 */
class PoiseuilleSolution2D : public Function<2> {
private:
  double H;    // Altezza del canale
  double Umax; // Velocità massima al centro del canale (y = H/2)

public:
  PoiseuilleSolution2D(double channel_height, double max_center_velocity)
      : H(channel_height), Umax(max_center_velocity) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    double y_norm = p.y / H;
    double ux = 4.0 * Umax * y_norm * (1.0 - y_norm);
    return utils::Vector<double, 2>{ux, 0.0};
  }
};

/**
 * Soluzione analitica per il flusso di Hagen-Poiseuille in un condotto
 * cilindrico (3D), asse parallelo a x.
 * Profilo parabolico di rivoluzione: u_x(r) = Umax * (1 - r^2/R^2), con
 * r la distanza dall'asse del tubo; u_y = u_z = 0.
 *
 * Fuori dal condotto (r >= R) vale zero, non il prolungamento negativo
 * della parabola: cosi' i nodi solidi della parete, dove il solver lascia
 * u = 0, non inquinano l'errore calcolato su tutta la griglia da
 * ErrorEvaluator<3>::integrate_difference().
 *
 * R e' il raggio *effettivo* della parete. Con bounce-back halfway la
 * parete non sta sui nodi solidi ma a meta' strada fra l'ultimo nodo di
 * fluido e il primo nodo solido, quindi vale R = r_inner + 0.5 se
 * r_inner e' il raggio passato a CylindricalShell.
 */
class HagenPoiseuilleSolution3D : public Function<3> {
private:
  double R;      // Raggio effettivo del condotto
  double Umax;   // Velocita' sull'asse
  double cy, cz; // Posizione dell'asse nel piano della sezione

public:
  HagenPoiseuilleSolution3D(double pipe_radius, double max_axis_velocity,
                            double axis_y, double axis_z)
      : R(pipe_radius), Umax(max_axis_velocity), cy(axis_y), cz(axis_z) {}

  utils::Vector<double, 3> value(const types::Coordinate<3> &p) const override {
    const double dy = p.y - cy;
    const double dz = p.z - cz;
    const double r2 = dy * dy + dz * dz;
    const double R2 = R * R;

    if (r2 >= R2) {
      return utils::Vector<double, 3>{0.0, 0.0, 0.0};
    }

    const double ux = Umax * (1.0 - r2 / R2);
    return utils::Vector<double, 3>{ux, 0.0, 0.0};
  }
};

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP
