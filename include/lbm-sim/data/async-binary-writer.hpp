#ifndef __LBM_SIM_DATA_ASYNC_BINARY_WRITER
#define __LBM_SIM_DATA_ASYNC_BINARY_WRITER

#include "lbm-sim/data/data-listener.hpp"

#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace lbm {

// Implementazione concreta di IDataListener: scrive in modo asincrono
// (su un thread dedicato) i chunk di dati ricevuti su un file binario.
// E' l'unico punto del codice che sa come/dove i dati finiscono su disco;
// i produttori (solver, LBMSimulation) non devono conoscere il formato
// del file, solo accodare byte grezzi tramite enqueueData.
class AsyncBinaryWriter : public IDataListener {
public:
  AsyncBinaryWriter(const std::string &path) : stop_(false) {
    file_.open(path, std::ios::out | std::ios::binary);
    if (!file_.is_open()) {
      std::cerr << "[ERROR] Impossibile aprire il file per la scrittura: "
                << path << std::endl;
    }
    worker_ = std::thread(&AsyncBinaryWriter::run, this);
  }

  ~AsyncBinaryWriter() override {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      stop_ = true;
    }
    cv_.notify_all();
    if (worker_.joinable()) {
      worker_.join();
    }
    if (file_.is_open()) {
      file_.flush();
      file_.close();
    }
  }

  // Implementazione di IDataListener::enqueueData.
  void acceptData(std::vector<char> data) override {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      queue_.push(std::move(data));
    }
    cv_.notify_one();
  }

private:
  void run() {
    while (true) {
      std::vector<char> chunk;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        // Attendi finche' non ci sono dati in coda o non arriva lo stop
        cv_.wait(lock, [this] { return stop_ || !queue_.empty(); });

        if (queue_.empty() && stop_) {
          break;
        }

        chunk = std::move(queue_.front());
        queue_.pop();
      }

      // Scrivi i dati su disco
      if (file_.is_open()) {
        file_.write(chunk.data(), chunk.size());
        file_.flush(); // Forziamo la scrittura fisica su disco ad ogni frame
      }
    }
  }

  std::ofstream file_;
  std::thread worker_;
  std::queue<std::vector<char>> queue_;
  std::mutex mtx_;
  std::condition_variable cv_;
  bool stop_;
};

} // namespace lbm

#endif // __LBM_SIM_DATA_ASYNC_BINARY_WRITER
