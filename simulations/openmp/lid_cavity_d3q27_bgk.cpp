// main_lidcavity_3d.cpp
//
// Lid-driven cavity 3D, collisione BGK, velocity set D3Q27, backend
// OpenMP.
//
// Compilazione (esempio):
//   g++ -std=c++17 -O3 -fopenmp -I<repo_root> main_lidcavity_3d.cpp \
//       -o lidcavity_3d
// (assumendo che <repo_root> contenga la cartella lbm-sim/, dato che
// tutti gli header sono inclusi come "lbm-sim/...").

#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"

#include "lbm/logging.hpp"

// QUILL LIB
#include "quill/LogMacros.h"

#include <memory>

namespace {

// Maschera dei confini per una cavita' cubica Nx*Ny*Nz: 5 pareti rigide
// (bounce-back semplice) + 1 parete mobile (il "lid", bounce-back con
// velocita' imposta) sulla faccia superiore z = Nz-1.
lbm::types::boundary_mask_t
build_boundary_mask(const lbm::types::DimPoint<3> &size) {
  lbm::types::boundary_mask_t mask(size.x * size.y * size.z, lbm::Solid::NONE);
  const lbm::Grid<3> grid(size);

#pragma omp parallel for collapse(3) schedule(static)
  for (std::size_t z = 0; z < size.z; ++z) {
    for (std::size_t y = 0; y < size.y; ++y) {
      for (std::size_t x = 0; x < size.x; ++x) {
        const bool on_boundary = (x == 0 || x == size.x - 1 || y == 0 ||
                                  y == size.y - 1 || z == 0 ||
                                  z == size.z - 1);
        if (!on_boundary)
          continue;

        const lbm::types::Coordinate<3> p(
            static_cast<int>(x), static_cast<int>(y), static_cast<int>(z));

        mask[grid.scalar_index(p)] = (z == size.z - 1)
                                         ? lbm::Solid::BB_MOVING_WALL
                                         : lbm::Solid::BB_RIGID_WALL;
      }
    }
  }

  return mask;
}

} // namespace

int main() {
  using namespace lbm;

  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  // --- Dominio: cubo Nx*Ny*Nz -------------------------------------
  // ATTENZIONE: con D3Q27 il costo per nodo e' 27/9 = 3x quello di
  // D2Q9, e il numero di nodi cresce come N^3 anziche' N^2: una
  // cavita' 128^3 ha gia' oltre 2M di nodi. Si consiglia di partire
  // piccoli (es. 32^3 o 48^3) per validare correttezza e prestazioni
  // prima di scalare.
  constexpr std::size_t N = 32;
  const types::DimPoint<3> grid_size(N, N, N);

  // --- Boundary mask: 5 pareti rigide + 1 lid mobile (z = N-1) -----
  const types::boundary_mask_t boundary_mask = build_boundary_mask(grid_size);

  // --- Parametri fisici ---------------------------------------------
  // init_vel: velocita' della parete mobile, diretta lungo x (dy=dz=0,
  // come da convenzione usata in Params per calcolare nu).
  const utils::Vector<double, 3> lid_velocity(0.1, 0.0, 0.0);
  const double reynolds = 100.0;

  const CollisionParams<3, CollisionModel::BGK> params(reynolds, grid_size,
                                                       lid_velocity);

  // --- Simulazione e solver -----------------------------------------
  constexpr unsigned int n_iters = 20000;
  constexpr unsigned int n_frames = 200; // ogni quante iter si salva

  LOG_INFO(
      main_logger,
      "Simulation #{} Parameters:\n\tGrid dimensions: {}\n\tReynolds number: "
      "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber of "
      "frames: {}",
      1, grid_size, reynolds, lid_velocity, n_iters, n_frames);

  LBMSimulation<3, D3Q27, CollisionModel::BGK> sim(grid_size, boundary_mask,
                                                   params);

  OpenMPSolver<3, D3Q27, CollisionModel::BGK> solver(n_iters, n_frames);

  // --- Output asincrono in .vtk (una serie temporale ParaView) -------
  // Lo stesso listener e' agganciato sia a `sim` (scrive l'header con
  // le dimensioni della griglia) sia a `solver` (scrive i frame delle
  // norme di velocita'): write_header() e write_norms() notificano due
  // DataObservable distinti (LBMSimulation e SolverBase), quindi vanno
  // attaccati entrambi allo stesso writer. VtkWriter usa il primo chunk
  // ricevuto (l'header) per conoscere nx,ny,nz e i successivi per
  // scrivere un file .vtk per frame in out/lidcavity3d_vtk/, piu' un
  // lidcavity3d.pvd che li raccoglie come serie temporale (apri quello
  // in ParaView, non i singoli .vtk).
  auto writer =
      std::make_shared<VtkWriter>("out/lidcavity3d_vtk", "lidcavity3d");
  sim.attachListener(writer);
  solver.attachListener(writer);

  sim.solve(solver);

  // NOTA: LBMSimulation::output() (estrazione della centerline)
  // e' ancora scritto assumendo un dominio
  // 2D (indicizza con size.x*y+x, ignora z) e non va usato cosi' com'e'
  // in 3D: per ora l'unico output affidabile in 3D e' la serie .vtk/.pvd
  // sopra. Va scritta una variante 3D (es. estrazione lungo z a x,y
  // fissati al centro) se serve un confronto quantitativo.

  return 0;
}
