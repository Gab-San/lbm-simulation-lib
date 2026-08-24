#ifndef __LBM_SIM_DATA_DATA_OBSERVABLE
#define __LBM_SIM_DATA_DATA_OBSERVABLE

#include "lbm-sim/data/data-listener.hpp"

#include <algorithm>
#include <memory>
#include <vector>

namespace lbm {

// Mixin per il ruolo di "observer/subject" nel pattern Observer/Listener.
// Solver e LBMSimulation ereditano (anche multiplo) da questa classe per
// poter notificare uno o piu' IDataListener (es. AsyncBinaryWriter) ogni
// volta che c'e' un nuovo chunk di dati da scrivere (header, norme, ...).
//
// Nota: non gestisce l'ownership dei listener, solo puntatori non
// proprietari — chi crea il listener (tipicamente LBMSimulation) resta
// responsabile della sua vita.

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
