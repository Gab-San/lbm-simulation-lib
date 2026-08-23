// main_lidcavity_3d.cpp
//
// Lid-driven cavity 3D, collisione BGK, velocity set D3Q27, backend
// OpenMP. Speculare al main 2D esistente (che usa LidCavity2D + D2Q9 +
// MPISolver2D), ma con la controparte 3D introdotta in questo commit:
// LidCavity3D (problems/problem_3d.hpp), D3Q27 (core/velocity-sets.hpp)
// e MPISolver3D (solver/solver-3d.hpp).
//
// Compilazione (esempio):
//   g++ -std=c++17 -O3 -fopenmp -I<repo_root> main_lidcavity_3d.cpp \
//       -o lidcavity_3d
// (assumendo che <repo_root> contenga la cartella lbm-sim/, dato che
// tutti gli header sono inclusi come "lbm-sim/...").

#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/problems/problem_3d.hpp"
#include "lbm-sim/solver/solver-3d.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

#include <memory>

int main() {
  using namespace lbm;

  // --- Dominio: cubo Nx*Ny*Nz -------------------------------------
  // ATTENZIONE: con D3Q27 il costo per nodo e' 27/9 = 3x quello di
  // D2Q9, e il numero di nodi cresce come N^3 anziche' N^2: una
  // cavita' 128^3 ha gia' oltre 2M di nodi. Si consiglia di partire
  // piccoli (es. 32^3 o 48^3) per validare correttezza e prestazioni
  // prima di scalare.
  constexpr std::size_t N = 32;
  const types::DimPoint<3> grid_size(N, N, N);

  // --- Boundary mask: 5 pareti rigide + 1 lid mobile (z = N-1) -----
  const types::boundary_mask_t boundary_mask =
      LidCavity3D::build_boundary_mask(grid_size);

  // --- Parametri fisici ---------------------------------------------
  // init_vel: velocita' della parete mobile, diretta lungo x (dy=dz=0,
  // come da convenzione usata in Params per calcolare nu).
  const utils::Vector<double, 3> lid_velocity(0.1, 0.0, 0.0);
  const double reynolds = 100.0;

  const Params<3, CollisionModel::BGK> params(reynolds, grid_size,
                                              lid_velocity);

  // --- Simulazione e solver -----------------------------------------
  constexpr unsigned int n_iters = 20000;
  constexpr unsigned int n_frames = 200; // ogni quante iter si salva

  LBMSimulation<3, D3Q27, CollisionModel::BGK> sim(grid_size, boundary_mask,
                                                    params);

  MPISolver3D<CollisionModel::BGK> solver(n_iters, n_frames);

  // --- Output asincrono su file binario ------------------------------
  // Lo stesso listener e' agganciato sia a `sim` (scrive l'header con
  // le dimensioni della griglia) sia a `solver` (scrive i frame delle
  // norme di velocita'): write_header() e write_norms() notificano due
  // DataObservable distinti (LBMSimulation e SolverBase), quindi vanno
  // attaccati entrambi allo stesso writer per finire nello stesso file.
  auto writer =
      std::make_shared<AsyncBinaryWriter>("out/lidcavity3d.bin");
  sim.attachListener(writer);
  solver.attachListener(writer);

  sim.solve(solver, LidCavity3D{});

  // NOTA: LBMSimulation::output() (estrazione della centerline per il
  // confronto con Ghia et al.) e' ancora scritto assumendo un dominio
  // 2D (indicizza con size.x*y+x, ignora z) e non va usato cosi' com'e'
  // in 3D: per ora l'unico output affidabile in 3D e' lo stream binario
  // via AsyncBinaryWriter sopra. Va scritta una variante 3D (es.
  // estrazione lungo z a x,y fissati al centro) se serve un confronto
  // quantitativo.

  return 0;
}
