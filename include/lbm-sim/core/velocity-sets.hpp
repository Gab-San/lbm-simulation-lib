#ifndef __LBM_SIM_CORE_VELOCITY_SETS_HPP
#define __LBM_SIM_CORE_VELOCITY_SETS_HPP

#include "lbm-sim/core/vector.hpp"

#ifdef __CUDACC__
#include "lbm-sim/backend/cuda-annotations.hpp"
#include "lbm-sim/backend/cuda-utils.cuh"

#include <cuda_runtime.h>
#endif

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {

namespace types {
template <unsigned short int dim> using VectorInt = utils::Vector<int, dim>;
}

struct D2Q9 {

  /// Number of dimensions
  static inline constexpr unsigned short int dim = 2;

  /// Number of directions
  static inline constexpr std::size_t ndir = 9;

  /// Weight in (dx,dy)=(0,0)
  static inline constexpr double w0 = 4.0 / 9.0;

  /// Weight for adjacent points
  static inline constexpr double ws = 1.0 / 9.0;

  /// Diagonal weight
  static inline constexpr double wd = 1.0 / 36.0;

  /**
   * Direction weights map.
   *
   * The direction numbering scheme is: \n
   * ------ + x \n
   * |7 4 8 \n
   * |3 0 1 \n
   * |6 2 5 \n
   * +\n
   * y
   */
  static inline constexpr double wi[ndir] = {w0, ws, ws, ws, ws,
                                             wd, wd, wd, wd};

  /**
   * Array of directions following the numbering scheme.
   *
   * The direction numbering scheme is: \n
   * ------ + x \n
   * |7 4 8 \n
   * |3 0 1 \n
   * |6 2 5 \n
   * \n
   * y
   */
  static inline constexpr types::VectorInt<2> dir[ndir] = {
      types::VectorInt<2>(0, 0),  types::VectorInt<2>(1, 0),
      types::VectorInt<2>(0, 1),  types::VectorInt<2>(-1, 0),
      types::VectorInt<2>(0, -1), types::VectorInt<2>(1, 1),
      types::VectorInt<2>(-1, 1), types::VectorInt<2>(-1, -1),
      types::VectorInt<2>(1, -1)};

  static inline constexpr std::size_t opp[ndir] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
};

// TODO: Uncomment when the time comes
//
// struct D3Q27 {
//
//   static constexpr unsigned short int dim = 3;
//   /// Number of directions
//   static constexpr std::size_t ndir = 27;
//   /// Weight in (dx,dy)=(0,0)
//   static constexpr double w0 = 8.0 / 27.0;
//   /// Weight for face points
//   static constexpr double wf = 2.0 / 27.0;
//   /// Weight for edge points
//   static constexpr double we = 1.0 / 54.0;
//   /// Weight for corner points
//   static constexpr double wc = 1.0 / 216.0;
//
//   // implemented in 3D but not used yet of directions and their opposite in
//   3D,
//   // following the numbering scheme:(x,y,z) 1(direction)
//   // --->2(opposite),3(direction) --->4(opposite),5(direction)
//   // --->6(opposite),7(direction) --->8(opposite),9(direction)
//   // --->10(opposite),11(direction) --->12(opposite),13(direction)
//   // --->14(opposite),15(direction) --->16(opposite),17(direction)
//   // --->18(opposite),19(direction) --->20(opposite),21(direction)
//   // --->22(opposite),23(direction) --->24(opposite),25(direction)
//   // --->26(opposite)
//   static LBM_CONST std::array<int, ndir> dirx = {
//       0,  1, -1, 0, 0, 0, 0,  1, -1, 1, -1, 1, -1, 1,
//       -1, 0, 0,  0, 0, 1, -1, 1, -1, 1, -1, 1, -1};
//
//   static LBM_CONST std::array<int, ndir> diry = {
//       0, 0, 0,  1, -1, 0, 0,  1, -1, -1, 1, 0,  0, 0,
//       0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1};
//
//   static LBM_CONST std::array<int, ndir> dirz = {
//       0, 0, 0,  0,  0, 1, -1, 0,  0, 0, 0,  1,  -1, -1,
//       1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1};
//
//   // FIXME: Check opposite
//   static LBM_CONST std::array<int, ndir> opp = {
//       0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
//       13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
// };
//

#ifdef __CUDACC__

namespace cuda {

template <unsigned short int dim, typename VelocitySet>
__constant__ types::VectorInt<dim> vs_dir[VelocitySet::ndir];

template <typename VelocitySet> __constant__ double vs_wi[VelocitySet::ndir];

template <typename VelocitySet>
__constant__ std::size_t
    vs_opp[VelocitySet::ndir]; // indice della direzione
                               // opposta (per il bounce-back)

template <unsigned short int dim, typename VelocitySet>
inline void upload_lattice_constants() {

  LBM_CUDA_CHECK(cudaMemcpyToSymbol(vs_dir<dim, VelocitySet>, VelocitySet::dir,
                                    sizeof(VelocitySet::dir)));
  LBM_CUDA_CHECK(cudaMemcpyToSymbol(vs_wi<VelocitySet>, VelocitySet::wi,
                                    sizeof(VelocitySet::wi)));
  LBM_CUDA_CHECK(cudaMemcpyToSymbol(vs_opp<VelocitySet>, VelocitySet::opp,
                                    sizeof(VelocitySet::opp)));
}

} // namespace cuda

#endif

} // namespace lbm

#endif // __LBM_SIM_CORE_VELOCITY_SETS_HPP
