#pragma once

#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/types/common.hpp"

#ifdef __CUDACC__
#include "lbm-sim/cuda/utils.cuh"

#include <cuda_runtime.h>
#endif

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {

struct D2Q9 {

  /// Number of dimensions
  static inline constexpr types::dim_t dim = 2;

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
   * |6 dim 5 \n
   * \n
   * y
   */
  static inline constexpr types::Direction<dim> dir[ndir] = {
      types::Direction<dim>(0, 0),  types::Direction<dim>(1, 0),
      types::Direction<dim>(0, 1),  types::Direction<dim>(-1, 0),
      types::Direction<dim>(0, -1), types::Direction<dim>(1, 1),
      types::Direction<dim>(-1, 1), types::Direction<dim>(-1, -1),
      types::Direction<dim>(1, -1)};

  static inline constexpr std::size_t opp[ndir] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
};

struct D3Q19 {

  static inline constexpr types::dim_t dim = 3;

  /// Number of directions
  static inline constexpr std::size_t ndir = 19;

  /// Weight in (dx,dy)=(0,0)
  static inline constexpr double w0 = 1 / 3.0;

  /// Weight for face points
  static inline constexpr double wf = 1 / 18.0;

  /// Weight for edge points
  static inline constexpr double we = 1 / 36.0;

  static inline constexpr double wi[ndir] = {
      w0 /*0*/, wf, wf, wf, wf, wf, wf /*6*/, we, we, we,
      we,       we, we, we, we, we, we,       we, we /*18*/};

  /**
   * Array of directions for D3Q19.
   *
   * \htmlonly
   * <iframe src="lbm_d3q19_directions.html" width="580" height="560"
   * style="border:none;"></iframe>
   * \endhtmlonly
   */
  static inline constexpr types::Direction<dim> dir[ndir] = {
      types::Direction<dim>(0, 0, 0) /*0*/,
      types::Direction<dim>(1, 0, 0),
      types::Direction<dim>(-1, 0, 0),
      types::Direction<dim>(0, 1, 0) /*3*/,
      types::Direction<dim>(0, -1, 0),
      types::Direction<dim>(0, 0, 1),
      types::Direction<dim>(0, 0, -1) /*6*/,
      types::Direction<dim>(1, 1, 0),
      types::Direction<dim>(-1, -1, 0),
      types::Direction<dim>(1, 0, 1),
      types::Direction<dim>(-1, 0, -1) /*10*/,
      types::Direction<dim>(0, 1, 1),
      types::Direction<dim>(0, -1, -1),
      types::Direction<dim>(1, -1, 0),
      types::Direction<dim>(-1, 1, 0) /*14*/,
      types::Direction<dim>(1, 0, -1),
      types::Direction<dim>(-1, 0, 1),
      types::Direction<dim>(0, 1, -1),
      types::Direction<dim>(0, -1, 1) /*18*/
  };

  static inline constexpr std::size_t opp[ndir] = {
      0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
};

struct D3Q27 {
  static inline constexpr types::dim_t dim = 3;

  /// Number of directions
  static constexpr std::size_t ndir = 27;
  /// Weight in (dx,dy)=(0,0)
  static constexpr double w0 = 8.0 / 27.0;
  /// Weight for face points
  static constexpr double wf = 2.0 / 27.0;
  /// Weight for edge points
  static constexpr double we = 1.0 / 54.0;
  /// Weight for corner points
  static constexpr double wc = 1.0 / 216.0;

  /**
   * Direction weights map.
   */
  static inline constexpr double wi[ndir] = {w0 /*0*/,  wf, wf, wf, wf, wf,
                                             wf /*6*/,  we, we, we, we, we,
                                             we,        we, we, we, we, we,
                                             we /*18*/, wc, wc, wc, wc, wc,
                                             wc,        wc, wc /*26*/};

  // FIXME: FIX DOCUMENTATION
  // implemented in 3D but not used yet of directions and their opposite in 3D,
  // following the numbering scheme:(x,y,z) 1(direction)
  // --->2(opposite),3(direction) --->4(opposite),5(direction)
  // --->6(opposite),7(direction) --->8(opposite),9(direction)
  // --->10(opposite),11(direction) --->12(opposite),13(direction)
  // --->14(opposite),15(direction) --->16(opposite),17(direction)
  // --->18(opposite),19(direction) --->20(opposite),21(direction)
  // --->22(opposite),23(direction) --->24(opposite),25(direction)
  // --->26(opposite)
  /**
   * Array of directions for D3Q27.
   *
   * \htmlonly
   * <iframe src="lbm_d3q27_directions.html" width="580" height="560"
   * style="border:none;"></iframe>
   * \endhtmlonly
   */
  static inline constexpr types::Direction<dim> dir[ndir] = {
      types::Direction<dim>(0, 0, 0) /*0*/,
      types::Direction<dim>(1, 0, 0) /*1*/,
      types::Direction<dim>(-1, 0, 0) /*2*/,
      types::Direction<dim>(0, 1, 0) /*3*/,
      types::Direction<dim>(0, -1, 0) /*4*/,
      types::Direction<dim>(0, 0, 1) /*5*/,
      types::Direction<dim>(0, 0, -1) /*6*/,
      types::Direction<dim>(1, 1, 0) /*7*/,
      types::Direction<dim>(-1, -1, 0) /*8*/,
      types::Direction<dim>(1, 0, 1) /*9*/,
      types::Direction<dim>(-1, 0, -1) /*10*/,
      types::Direction<dim>(0, 1, 1) /*11*/,
      types::Direction<dim>(0, -1, -1) /*12*/,
      types::Direction<dim>(1, -1, 0) /*13*/,
      types::Direction<dim>(-1, 1, 0) /*14*/,
      types::Direction<dim>(1, 0, -1) /*15*/,
      types::Direction<dim>(-1, 0, 1) /*16*/,
      types::Direction<dim>(0, 1, -1) /*17*/,
      types::Direction<dim>(0, -1, 1) /*18*/,
      types::Direction<dim>(1, 1, 1) /*19*/,
      types::Direction<dim>(-1, -1, -1) /*20*/,
      types::Direction<dim>(1, 1, -1) /*21*/,
      types::Direction<dim>(-1, -1, 1) /*22*/,
      types::Direction<dim>(1, -1, 1) /*23*/,
      types::Direction<dim>(-1, 1, -1) /*24*/,
      types::Direction<dim>(-1, 1, 1) /*25*/,
      types::Direction<dim>(1, -1, -1) /*26*/};

  static inline constexpr std::size_t opp[ndir] = {
      0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
      13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
};

#ifdef __CUDACC__

namespace cuda {

template <unsigned short int dim, typename VelocitySet>
__constant__ types::Direction<dim> vs_dir[VelocitySet::ndir];

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
