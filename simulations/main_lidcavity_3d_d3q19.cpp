// main_lidcavity_3d_d3q19.cpp
//
// Identico a main_lidcavity_3d.cpp, ma con il velocity set D3Q19 al
// posto di D3Q27 (19 direzioni invece di 27: niente vicini di
// vertice/corner, solo riposo + facce + spigoli). Costa meno per nodo
// (19 vs 27 popolazioni da streammare/collidere) a scapito di un po'
// di isotropia rispetto a D3Q27 - per la lid cavity standard e' comunque
// una scelta molto comune in letteratura.
//
// MPISolver3D e SolverBase3D sono ora templati su VelocitySet, quindi
// l'unica differenza col main D3Q27 e' il tipo passato a
// LBMSimulation<...> e MPISolver3D<...>.
//
// Compilazione: identica a main_lidcavity_3d.cpp.

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

  LBMSimulation<3, D3Q19, CollisionModel::BGK> sim(grid_size, boundary_mask,
                                                    params);

  MPISolver3D<CollisionModel::BGK, D3Q19> solver(n_iters, n_frames);

  // --- Output asincrono su file binario ------------------------------
  // Lo stesso listener e' agganciato sia a `sim` (scrive l'header con
  // le dimensioni della griglia) sia a `solver` (scrive i frame delle
  // norme di velocita'): write_header() e write_norms() notificano due
  // DataObservable distinti (LBMSimulation e SolverBase), quindi vanno
  // attaccati entrambi allo stesso writer per finire nello stesso file.
  auto writer =
      std::make_shared<AsyncBinaryWriter>("out/lidcavity3d_d3q19.bin");
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
