#ifndef __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP
#define __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP

#include "lbm-sim/analysis/error.hpp"

#include <memory>
#include <stdexcept>

namespace lbm {
namespace analysis {

/**
 * Soluzione analitica per il Flusso di Couette (2D).
 * Parete inferiore fissa (y = 0), parete superiore mobile (y = H) con velocità Umax.
 * Profilo di velocità lineare: u_x(y) = Umax * (y / H), u_y = 0.
 */
class CouetteSolution2D : public Function<2> {
private:
  double H;    // Altezza del canale
  double Umax; // Velocità della parete mobile superiore

public:
  CouetteSolution2D(double channel_height, double max_velocity)
      : H(channel_height), Umax(max_velocity) {}

  utils::Vector<double, 2>
  value(const types::ValuePoint<2> &p) const override {
    // p.y rappresenta la coordinata verticale[cite: 3]
    double ux = Umax * (p.y / H);
    return utils::Vector<double, 2>{ux, 0.0};
  }
};

/**
 * Soluzione analitica per il Flusso di Poiseuille (2D).
 * Guidato da gradiente di pressione / forza di volume tra due pareti fisse (y = 0 e y = H).
 * Profilo di velocità parabolico: u_x(y) = 4 * Umax * (y / H) * (1 - y / H), u_y = 0.
 */
class PoiseuilleSolution2D : public Function<2> {
private:
  double H;    // Altezza del canale
  double Umax; // Velocità massima al centro del canale (y = H/2)

public:
  PoiseuilleSolution2D(double channel_height, double max_center_velocity)
      : H(channel_height), Umax(max_center_velocity) {}

  utils::Vector<double, 2>
  value(const types::ValuePoint<2> &p) const override {
    double y_norm = p.y / H;
    double ux = 4.0 * Umax * y_norm * (1.0 - y_norm);
    return utils::Vector<double, 2>{ux, 0.0};
  }
};

enum class FlowType { Couette, Poiseuille };

/**
 * Factory: costruisce la Function<2> giusta in base a flow_type,
 * riusando H/Umax presi dalla stessa Config usata per lanciare la
 * simulazione (channel_height tipicamente grid_size.y - 1, ref_velocity
 * tipicamente init_vel.dx). Non normalizza nulla: passa esattamente i
 * valori fisici che vuoi confrontare.
 */
inline std::unique_ptr<Function<2>>
make_exact_solution(FlowType flow_type, double channel_height,
                    double ref_velocity) {
  switch (flow_type) {
  case FlowType::Couette:
    return std::make_unique<CouetteSolution2D>(channel_height, ref_velocity);
  case FlowType::Poiseuille:
    return std::make_unique<PoiseuilleSolution2D>(channel_height,
                                                   ref_velocity);
  }
  throw std::invalid_argument(
      "make_exact_solution(): unhandled FlowType value");
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP