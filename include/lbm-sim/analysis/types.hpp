#ifndef __LBM_SIM_ANALYSIS_TYPES_HPP
#define __LBM_SIM_ANALYSIS_TYPES_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {
namespace analysis {

template <types::dim_t dim> class Function {
public:
  virtual ~Function() = default;

  // Restituisce la velocità vettoriale esatta nel punto continuo P
  virtual utils::Vector<double, dim>
  value(const types::Coordinate<dim> &p) const = 0;
};

enum class NormType { L1, L2, L2_squared, Linfty };

inline const char *to_string(NormType t) {
  switch (t) {
  case NormType::L1:
    return "L1";
  case NormType::L2:
    return "L2";
  case NormType::L2_squared:
    return "L2^2";
  case NormType::Linfty:
    return "Linf";
  }
  return "unknown";
}

/**
 * Risultato di un confronto con una soluzione analitica generica
 * (Couette, Poiseuille, ...) dove la norma è scelta dal chiamante.
 */
struct NormErrorResult {
  double relative;
  double absolute;
  NormType norm_type;

  // Metriche aggiuntive per una lettura piu' immediata dell'errore.
  // I campi normalized sono rapporti rispetto alla velocita' di riferimento
  // del problema (per Ghia i dati sono gia' normalizzati con U_lid).
  double rmse = 0.0;
  double linf = 0.0;
  double rmse_normalized = 0.0;
  double linf_normalized = 0.0;
  std::size_t sample_count = 0;
};

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_TYPES_HPP
