/**
 * @file data-listener.hpp
 * @brief IDataListener: the one-method interface every output sink
 *        implements.
 *
 * The narrowest possible contract -- a chunk of bytes, nothing else -- so
 * that the producers (LBMSimulation, the solvers) know nothing about file
 * formats, and a new sink is one class with one method.
 *
 * The order and meaning of the chunks is documented on the "Output formats"
 * page: a grid header first, then one frame per emission.
 */

#ifndef __LBM_SIM_DATA_DATA_LISTENER
#define __LBM_SIM_DATA_DATA_LISTENER

#include <vector>

namespace lbm {

/**
 * @brief An object that recieves data to be handled and outputted
 * is denoted as a _Data Listener_.
 */
class IDataListener {
public:
  virtual ~IDataListener() = default;

  /**
   * Recieves a vector of data to be handled.
   *
   * @param[in] data  chunk of data to be handled
   */
  virtual void acceptData(std::vector<char> data) = 0;
};

} // namespace lbm

#endif // __LBM_SIM_DATA_DATA_LISTENER
