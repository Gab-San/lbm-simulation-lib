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

struct D3Q27 {
    
    static constexpr int dim = 3;
    /// Number of directions
    static constexpr std::size_t ndir = 27;
    /// Weight in (dx,dy)=(0,0)
    static constexpr double w0 = 8.0/27.0;
    /// Weight for face points
    static constexpr double wf = 2.0/27.0;
    /// Weight for edge points
    static constexpr double we = 1.0/54.0;
    /// Weight for corner points
    static constexpr double wc = 1.0/216.0;
 
 
    // implemented in 3D but not used yet of directions and their opposite in 3D, following the numbering scheme:(x,y,z)
    // 1(direction) --->2(opposite),3(direction) --->4(opposite),5(direction) --->6(opposite),7(direction) --->8(opposite),9(direction) --->10(opposite),11(direction) --->12(opposite),13(direction) --->14(opposite),15(direction) --->16(opposite),17(direction) --->18(opposite),19(direction) --->20(opposite),21(direction) --->22(opposite),23(direction) --->24(opposite),25(direction) --->26(opposite)
    static constexpr std::array<int, ndir3D> dirx = {
        0,
        1,-1, 0, 0, 0, 0,
        1,-1, 1,-1, 1,-1, 1,-1,
        0, 0, 0, 0,
        1,-1, 1,-1, 1,-1, 1,-1
    };

    static constexpr std::array<int, ndir3D> diry = {
        0,
        0, 0, 1,-1, 0, 0,
        1,-1,-1, 1, 0, 0, 0, 0,
        1,-1, 1,-1,
        1,-1, 1,-1,-1, 1,-1, 1
    };

    static constexpr std::array<int, ndir3D> dirz = {
        0,
        0, 0, 0, 0, 1,-1,
        0, 0, 0, 0, 1,-1,-1, 1,
        1,-1,-1, 1,
        1,-1,-1, 1, 1,-1,-1, 1
    };
    
    static constexpr std::array<int, ndir3D> opp = {
        0,
        2, 1, 4, 3, 6, 5,
        8, 7,10, 9,12,11,14,13,
        16,15,18,17,
        20,19,22,21,24,23,26,25
    };  

}