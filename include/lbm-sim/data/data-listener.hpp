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
