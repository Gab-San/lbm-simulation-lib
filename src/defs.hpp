#ifndef DEFS_HPP
#define DEFS_HPP

/**
 * \file
 * This file contains some definitions used
 * by all the code implementing the lid driven cavity
 * simulation.
 */
namespace lbm_lbm {

    /**
    * This enumeration serves to identify the type of 
    * Collision Operator used in the LBM algorithm.
    *
    * \important TRT is NOT currently supported
    */
    enum CollisionOperator {
	/** Two-Relaxation-Times */TRT, 
	/** Bhatnagar-Gross-Krook */BGK
    };
}

#endif
