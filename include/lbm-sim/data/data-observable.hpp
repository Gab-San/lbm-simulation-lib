/**
 * @file data-observable.hpp
 * @brief DataObservable: the subject half of the listener pattern.
 *
 * Mixed into LBMSimulation and SolverBase, which are two *separate*
 * observables: the simulation emits the grid header, the solver emits the
 * frames. A listener that needs the whole stream has to be attached to both.
 */

#ifndef __LBM_SIM_DATA_DATA_OBSERVABLE
#define __LBM_SIM_DATA_DATA_OBSERVABLE

#include "lbm-sim/data/data-listener.hpp"

#include <algorithm>
#include <memory>
#include <vector>

namespace lbm {

// Mixin providing the "subject" role of the Observer/Listener pattern.
// Solver and LBMSimulation inherit from this class (possibly through
// multiple inheritance) so they can notify one or more IDataListener
// (e.g. AsyncBinaryWriter) whenever a new chunk of data has to be written
// (header, norms, ...).
//
// Note: it does not own the listeners, it only keeps non-owning handles --
// whoever creates the listener (typically LBMSimulation) stays responsible
// for its lifetime.

/// Attach/detach and notify are NOT thread-safe with respect to each other.
// Contract: all attachListener()/detachListener() calls must happen
// before dispatch begins and after it has fully stopped for a given
// simulation run — never concurrently with notifyListeners().
class DataObservable {
public:
  virtual ~DataObservable() = default;

  /**
   * @brief Attaches a listener
   *
   * After calling this function, everytime the observable dispatches
   * an event the attached listener will be notified of it.
   *
   * @param[in] listener The listener to be attached
   */
  void attachListener(std::shared_ptr<IDataListener> listener) {
    if (listener != nullptr) {
      listeners_.push_back(listener);
    }
  }

  /**
   * @brief Detaches a listener
   *
   * After calling this function, dispached events will be not be
   * handed to the detached listener.
   *
   * @param listener The listener to detach
   */
  void detachListener(std::shared_ptr<IDataListener> listener) {
    listeners_.erase(
        std::remove(listeners_.begin(), listeners_.end(), listener),
        listeners_.end());
  }

protected:
  /**
   * @brief Notifies attached listeners of an event dispatch
   *
   * The data buffer is copied for all listeners, except for the last,
   * for which it is moved.
   *
   * @param data The data to be dispatched
   */
  void notifyListeners(std::vector<char> &&data) const {
    for (std::size_t i = 0; i < listeners_.size(); ++i) {
      if (i + 1 == listeners_.size()) {
        listeners_[i]->acceptData(std::move(data));
      } else {
        listeners_[i]->acceptData(data);
      }
    }
  }

private:
  std::vector<std::shared_ptr<IDataListener>> listeners_;
};

} // namespace lbm

#endif // __LBM_SIM_DATA_DATA_OBSERVABLE
