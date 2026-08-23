#ifndef __LBM_SIM_SOLVER_SOLVER_3D
#define __LBM_SIM_SOLVER_SOLVER_3D

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/backend/metadata.hpp"
#include "lbm-sim/boundaries.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"

#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

// C++ STANDARD LIB
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

// OMP LIB
#include <omp.h>

namespace lbm {

// Speculare a MPISolver2D (solver/solver-2d.hpp): stessa struttura
// stream+collide, stesso switch sulle boundary condition, ma con un
// terzo indice di griglia (z). Templata su VelocitySet (default D3Q27)
// cosi' lo stesso codice serve sia per D3Q27 sia per D3Q19 - la logica
// di streaming/collisione non dipende dal numero di direzioni, solo
// VelocitySet::ndir/dir/wi/opp cambiano. Se in solver-2d.hpp cambia la
// logica di streaming/collisione, la stessa modifica va replicata qui.
template <enum CollisionModel cm_t, typename VelocitySet = D3Q27>
class MPISolver3D : public SolverBase3D<VelocitySet, cm_t, ExecutionBackend::OPEN_MP> {
  static_assert(VelocitySet::dim == 3,
               "MPISolver3D richiede un VelocitySet con dim == 3 "
               "(es. D3Q27, D3Q19)");
  using Base = SolverBase3D<VelocitySet, cm_t, ExecutionBackend::OPEN_MP>;
  static constexpr unsigned short int DIM = 3;

public:
  MPISolver3D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {};

  virtual ~MPISolver3D() = default;

  void init_equilibrium(const Lattice<3> &lattice,
                        std::vector<double> &part_stream) const override {
    using utils::ops::dot;
#pragma omp parallel for shared(lattice, part_stream) schedule(static)         \
    collapse(3)
    for (unsigned int z = 0; z < lattice.grid.size.z; ++z) {
      for (unsigned int y = 0; y < lattice.grid.size.y; ++y) {
        for (unsigned int x = 0; x < lattice.grid.size.x; ++x) {
          const types::Coordinate<3> p(x, y, z);

          double r = lattice.rho[lattice.grid.scalar_index(p)];
          const utils::Vector<double, 3> u =
              lattice.u[lattice.grid.scalar_index(p)];
          const double u_sq = dot(u, u);

#pragma omp simd
          for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
            double cidotu = dot(VelocitySet::dir[i], u);

            part_stream[lattice.grid.field_index(p, i, VelocitySet::ndir)] =
                VelocitySet::wi[i] * r *
                (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
          }
        }
      }
    }
  };

  void solve(Lattice<3> &lattice, const Params<3, cm_t> &params_,
             std::vector<double> &ffrom,
             std::vector<double> &fto) const override {
    const CollisionStrategy<3, VelocitySet, cm_t, OPEN_MP> cs(params_);
    for (unsigned int iter = 0; iter < this->niters; iter++) {
      bool save = (iter % this->nskips == 0);
      update_stream_collide(lattice, cs, ffrom, fto, save);
      std::swap(ffrom, fto);
      if (save)
        write_norms(lattice);
    }
  }

private:
  void update_stream_collide(
      Lattice<3> &lattice, const CollisionStrategy<3, VelocitySet, cm_t, OPEN_MP> &cs,
      const std::vector<double> &ffrom, std::vector<double> &fto, bool save,
      const ExecutionContext<OPEN_MP> &context =
          ExecutionContext<OPEN_MP>{}) const {

#pragma omp parallel for shared(ffrom, fto, cs, lattice, save)                 \
    schedule(static) collapse(3)
    for (std::size_t z = 0; z < lattice.grid.size.z; ++z) {
      for (std::size_t y = 0; y < lattice.grid.size.y; ++y) {
        for (std::size_t x = 0; x < lattice.grid.size.x; ++x) {
          std::array<double, VelocitySet::ndir> fp;
          const utils::Point<int, 3> p(x, y, z);

          double r_wall = 0.0;
#pragma omp unroll full
          for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
            r_wall += ffrom[lattice.grid.field_index(p, i, VelocitySet::ndir)];
          }

#pragma omp unroll full
          for (unsigned int diridx = 0; diridx < VelocitySet::ndir; ++diridx) {
            const types::Coordinate<3> src = p - VelocitySet::dir[diridx];

            if (!lattice.grid.contains(src)) {
              apply_boundary_conditions(lattice.boundary_mask, fp, ffrom,
                                        diridx, lattice, p, r_wall,
                                        cs.params.init_vel, context);
            } else {
              fp[diridx] =
                  ffrom[lattice.grid.field_index(src, diridx, VelocitySet::ndir)];
            }
          }

          // rho = sum_i fi ; rho*u = sum_i fi * ci
          double r = 0.0;
          utils::Vector<double, 3> u(0, 0, 0);

#pragma omp unroll full
          for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
            r += fp[i];
            u += VelocitySet::dir[i] * fp[i];
          }

          u /= r;

          const auto s_idx = lattice.grid.scalar_index(p);
          if (save ||
              lattice.boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_INLET ||
              lattice.boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_OUTLET) {
            lattice.rho[s_idx] = r;
            lattice.u[s_idx] = u;
          }

          cs.apply(p, u, r, fp, lattice.grid);

#pragma omp simd
          for (auto i = 0; i < VelocitySet::ndir; i++) {
            fto[lattice.grid.field_index(p, i, VelocitySet::ndir)] = fp[i];
          }
        }
      }
    }
  };

  void apply_boundary_conditions(const types::boundary_mask_t &boundary_mask,
                                 std::array<double, VelocitySet::ndir> &fp,
                                 const std::vector<double> &ffrom,
                                 const std::size_t diridx,
                                 const Lattice<3> &lattice,
                                 const types::Coordinate<3> p,
                                 const double &localrho,
                                 const utils::Vector<double, 3> u0,
                                 const ExecutionContext<OPEN_MP> &context =
                                     ExecutionContext<OPEN_MP>{}) const {
    (void)context;

    types::boundary_t b = boundary_mask[lattice.grid.scalar_index(p)];

    switch (b) {
    case Solid::BB_RIGID_WALL:
      Solid::apply_bb_rigid_wall<3, VelocitySet>(fp, ffrom, diridx, lattice.grid, p);
      break;
    case Solid::BB_MOVING_WALL:
      Solid::apply_bb_moving_wall<3, VelocitySet>(fp, ffrom, diridx, lattice.grid, p,
                                            localrho, u0);
      break;
    case Solid::PERIODIC:
      Solid::apply_periodic<3, VelocitySet>(fp, ffrom, diridx, lattice.grid, p);
      break;
    case Solid::PRESSURE_PERIODIC_INLET:
      Solid::apply_periodic_with_pressure_variation<3, VelocitySet>(
          fp.data(), ffrom.data(), diridx, lattice.grid, p, lattice.rho.data(),
          lattice.pin, lattice.u.data());
      break;
    case Solid::PRESSURE_PERIODIC_OUTLET:
      Solid::apply_periodic_with_pressure_variation<3, VelocitySet>(
          fp.data(), ffrom.data(), diridx, lattice.grid, p, lattice.rho.data(),
          lattice.pout, lattice.u.data());
      break;
    default:
      break;
    }
  }

  void write_norms(const Lattice<3> &lattice) const override {
    using utils::ops::dot;
    std::vector<float> vsq(lattice.grid.getArea());

#pragma omp parallel for collapse(3)
    for (unsigned int z = 0; z < lattice.grid.size.z; ++z) {
      for (unsigned int y = 0; y < lattice.grid.size.y; ++y) {
        for (unsigned int x = 0; x < lattice.grid.size.x; ++x) {
          const types::Coordinate<3> p(x, y, z);
          const auto &vel = lattice.u[lattice.grid.scalar_index(p)];

          vsq[lattice.grid.scalar_index(p)] =
              static_cast<float>(std::sqrt(dot(vel, vel)));
        }
      }
    }

    std::vector<char> buf(vsq.size() * sizeof(float));
    std::memcpy(buf.data(), vsq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }

}; // class MPISolver3D

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_3D
