#pragma once

#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/base.hpp"
#include <array>
#include <stdexcept>
#include <string>

namespace lbm::config {

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
template <types::dim_t dim> struct SimulationConfig {

  /// [solver].backend: "openmp" o "cuda".
  ExecutionBackend backend = ExecutionBackend::OPEN_MP;

  /// Operatore di collisione.
  CollisionModel collision = CollisionModel::BGK;

  /// Celle della griglia. `nz` vale 1 (e non va letto) quando dim == 2.
  std::array<uint64_t, dim> grid_size;

  std::string name;
  /// Numero di Reynolds.
  double reynolds = 100.0;

  /// Velocita' di riferimento. `init_vel_z` esiste solo per dim == 3.
  ///
  /// \note init_vel_x non e' solo la velocita' iniziale: e' la velocita'
  /// caratteristica con cui CollisionParams calcola nu = u*Ny/Re, quindi
  /// deve essere strettamente positiva anche nei problemi (Poiseuille) in
  /// cui non muove nessuna parete.
  utils::Vector<double, dim> u0;

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
};

} // namespace lbm::config
