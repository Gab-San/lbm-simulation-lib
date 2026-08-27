#ifndef __LBM_SIM_CONFIG_SIMULATION_CONFIG_HPP
#define __LBM_SIM_CONFIG_SIMULATION_CONFIG_HPP

#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/base.hpp"
#include "lbm-sim/types/common.hpp"

#include <stdexcept>
#include <string>

namespace lbm {
namespace config {

/**
 * \brief Errore di configurazione.
 *
 * Copre tutto cio' che rende una configurazione inutilizzabile *prima* che
 * la simulazione parta: file illeggibile, TOML malformato, campo
 * obbligatorio mancante, valore fuori range, oppure una configurazione
 * coerente in se' ma incompatibile con il binario che la sta leggendo
 * (dimensione, tipo di problema, operatore di collisione, backend).
 *
 * E' un tipo a parte e non un std::runtime_error generico proprio perche'
 * i main possano distinguere "l'utente ha passato la config sbagliata"
 * (si stampa il messaggio e si esce con 1) da un errore che nasce durante
 * il solve.
 */
class ConfigError : public std::runtime_error {
public:
  explicit ConfigError(const std::string &what) : std::runtime_error(what) {}
};


/**
 * \brief Configurazione di *una* simulazione, letta da un singolo file
 * .toml piu' i due path di output presi dalla riga di comando.
 *
 * La dimensione del problema non e' un campo libero: la si deduce dalla
 * presenza di [grid].nz (assente -> 2D, presente -> 3D) e la si puo'
 * dichiarare esplicitamente con [grid].dim per farla ricontrollare dal
 * parser. Vedi config-parser.hpp.
 *
 * @verbatim
   [problem]
   type = "lid_cavity"

   [grid]
   nx = 129
   ny = 129
   # nz = 129     # presente solo per i problemi 3D

   [physics]
   collision   = "BGK"
   reynolds    = 100.0
   init_vel_x  = 0.1

   [solver]
   backend = "openmp"
   niters  = 10000
   nframes = 100
   @endverbatim
 */
struct SimulationConfig {

  /// Nome della configurazione: il nome del file .toml senza estensione.
  /// Serve solo per i messaggi di log e di errore.
  std::string name;

  /// [problem].type: "lid_cavity", "couette", "poiseuille", ...
  /// La geometria vera e propria (pareti, inlet/outlet) resta nei main:
  /// qui c'e' solo l'etichetta che dice a quale main appartiene la config.
  std::string problem;

  /// [solver].backend: "openmp" o "cuda".
  std::string backend = "openmp";

  /// Dimensione del problema, dedotta dal parser: 2 o 3.
  types::dim_t dim = 2;

  /// Celle della griglia. `nz` vale 1 (e non va letto) quando dim == 2.
  unsigned int nx = 0;
  unsigned int ny = 0;
  unsigned int nz = 1;

  /// Operatore di collisione.
  CollisionModel collision = CollisionModel::BGK;

  /// Numero di Reynolds.
  double reynolds = 100.0;

  /// Velocita' di riferimento. `init_vel_z` esiste solo per dim == 3.
  ///
  /// \note init_vel_x non e' solo la velocita' iniziale: e' la velocita'
  /// caratteristica con cui CollisionParams calcola nu = u*Ny/Re, quindi
  /// deve essere strettamente positiva anche nei problemi (Poiseuille) in
  /// cui non muove nessuna parete.
  double init_vel_x = 0.0;
  double init_vel_y = 0.0;
  double init_vel_z = 0.0;

  /// Numero di iterazioni del solver.
  unsigned int niters = 0;

  /// Numero di frame salvati durante il solve. Ogni frame contiene le
  /// norme della velocita' a un dato passo temporale.
  unsigned int nframes = 0;

  /// File di output delle norme (secondo argomento della riga di comando).
  std::string frames_out;

  /// File di output del profilo di velocita' (terzo argomento della riga
  /// di comando).
  std::string profile_out;

  /// Dimensioni della griglia nel tipo usato dalla libreria.
  /// `dim` va scelto a compile time perche' DimPoint<2> e DimPoint<3> sono
  /// tipi distinti; ensure_dim() garantisce che coincida con this->dim.
  template <types::dim_t dim_> types::DimPoint<dim_> grid_size() const {
    static_assert(dim_ == 2 || dim_ == 3,
                  "grid_size(): supportate solo 2 e 3 dimensioni");
    if constexpr (dim_ == 2) {
      return types::DimPoint<2>(nx, ny);
    } else {
      return types::DimPoint<3>(nx, ny, nz);
    }
  }

  /// Velocita' di riferimento nel tipo usato dalla libreria.
  template <types::dim_t dim_> utils::Vector<double, dim_> velocity() const {
    static_assert(dim_ == 2 || dim_ == 3,
                  "velocity(): supportate solo 2 e 3 dimensioni");
    if constexpr (dim_ == 2) {
      return utils::Vector<double, 2>(init_vel_x, init_vel_y);
    } else {
      return utils::Vector<double, 3>(init_vel_x, init_vel_y, init_vel_z);
    }
  }

  /// Numero totale di celle: nx*ny in 2D, nx*ny*nz in 3D.
  std::size_t total_cells() const {
    const std::size_t base =
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    return (dim == 3) ? base * static_cast<std::size_t>(nz) : base;
  }
};

} // namespace config
} // namespace lbm

#endif // __LBM_SIM_CONFIG_SIMULATION_CONFIG_HPP
