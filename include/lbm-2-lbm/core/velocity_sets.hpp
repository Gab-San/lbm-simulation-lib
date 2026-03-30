#pragma once

#include <cstddef>
#include <array>

namespace lbm {

struct D2Q9 {
    
    static constexpr int dim = 2;

    /// Number of directions
    static constexpr std::size_t ndir = 9;

    /// Weight in (dx,dy)=(0,0)
    static constexpr double w0 = 4.0/9.0;

    /// Weight for adjacent points
    static constexpr double ws = 1.0/9.0;

    /// Diagonal weight
    static constexpr double wd = 1.0/36.0;


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
    static constexpr std::array<double, ndir> wi = {w0, ws, ws, ws, ws, wd, wd, wd, wd};

    /**
     * Array of directions in the x direction following the numbering scheme.
     *
     * The direction numbering scheme is: \n
     * ------ + x \n
     * |7 4 8 \n
     * |3 0 1 \n
     * |6 2 5 \n
     * \n
     * y
     */
    static constexpr std::array<int, ndir> dirx = {0,1,0,-1,0,1,-1,-1,1};

    /**
     * Array of direction in the y direction following the numbering scheme.
     *
     * The direction numbering scheme is: \n
     * ------ + x \n
     * |7 4 8 \n
     * |3 0 1 \n
     * |6 2 5 \n
     * \n
     * y
     */
    static constexpr std::array<int, ndir> diry = {0,0,1,0,-1,1,1,-1,-1};

    static constexpr std::array<int, ndir> opp = {0, 3, 4, 1, 2, 6, 5, 7, 8};
};

}
