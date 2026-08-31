/**
 * @file async-binary-writer.hpp
 * @brief AsyncBinaryWriter: the default sink, a raw binary file written off
 *        the simulation thread.
 *
 * Chunks are queued and drained by a dedicated worker, so a slow disk stalls
 * the writer rather than the solver. The file is the concatenation of the
 * chunks exactly as they were handed over -- see the "Output formats" page
 * for the layout.
 *
 * @warning The destructor is what drains the queue and closes the file. A
 *          writer destroyed before the run ends loses the tail; one that
 *          leaks is never flushed at all.
 */

#ifndef __LBM_SIM_DATA_ASYNC_BINARY_WRITER
#define __LBM_SIM_DATA_ASYNC_BINARY_WRITER

#include "lbm-sim/data/data-listener.hpp"

#include "lbm-sim/logging.hpp"

#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace lbm {

// Concrete implementation of IDataListener: asynchronously writes the data
// chunks it receives to a binary file, on a dedicated thread.
// It is the only place in the code that knows how/where the data ends up on
// disk; the producers (solver, LBMSimulation) need not know the file format,
// only enqueue raw bytes through enqueueData.
class AsyncBinaryWriter : public IDataListener {
  const std::string path;

public:
  AsyncBinaryWriter(const std::string &path) : path(path), stop_(false) {
    logging::Logger *writer_logger = logging::create_or_get_logger("writer");
    std::filesystem::path frames_path_out(path);
    std::filesystem::path output_filename =
        frames_path_out.parent_path() / frames_path_out.stem();
    output_filename += ".bin";
    file_.open(output_filename, std::ios::out | std::ios::binary);
    if (!file_.is_open()) {
      LBM_LOG_CRITICAL(writer_logger, "Cannot open {} for writing", path);
      throw std::invalid_argument("Can't open " + output_filename.string() +
                                  " for writing");
    }
    worker_ = std::thread(&AsyncBinaryWriter::run, this);
    LBM_LOG_DEBUG(writer_logger, "File {} opened for binary writing", path);
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

  // Implementation of IDataListener::enqueueData.
  void acceptData(std::vector<char> data) override {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      queue_.push(std::move(data));
    }
    cv_.notify_one();
  }

private:
  void run() {
    size_t writes_since_flush = 0;
    constexpr size_t flush_interval = 32;

    while (true) {
      std::vector<char> chunk;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        // Wait until there is data in the queue or a stop is requested.
        cv_.wait(lock, [this] { return stop_ || !queue_.empty(); });

        if (queue_.empty() && stop_) {
          break;
        }

        chunk = std::move(queue_.front());
        queue_.pop();
      }

      // Write the data to disk.
      if (file_.is_open()) {
        file_.write(chunk.data(), chunk.size());

        if (++writes_since_flush >= flush_interval) {
          file_.flush();
          writes_since_flush = 0;
        }
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
