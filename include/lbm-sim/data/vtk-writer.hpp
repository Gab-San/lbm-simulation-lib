#ifndef __LBM_SIM_DATA_VTK_WRITER
#define __LBM_SIM_DATA_VTK_WRITER

#include "lbm-sim/data/data-listener.hpp"

#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace lbm {

class VtkWriter : public IDataListener {
public:
  VtkWriter(const std::string &output_dir, const std::string &basename)
      : output_dir_(output_dir), basename_(basename), have_dims_(false), nx_(0),
        ny_(0), nz_(1), frame_count_(0), stop_(false) {
    std::filesystem::create_directories(output_dir_);
    worker_ = std::thread(&VtkWriter::run, this);
  }

  explicit VtkWriter(const std::string &path)
      : VtkWriter(dir_from_path(path), basename_from_path(path)) {}

  ~VtkWriter() override {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      stop_ = true;
    }
    cv_.notify_all();
    if (worker_.joinable()) {
      worker_.join();
    }
  }

  void acceptData(std::vector<char> data) override {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      queue_.push(std::move(data));
    }
    cv_.notify_one();
  }

private:
  /// "out/norms_lid_cavity.bin" -> "norms_lid_cavity"
  static std::string basename_from_path(const std::string &path) {
    return std::filesystem::path(path).stem().string();
  }

  /// "out/norms_lid_cavity.bin" -> "out/norms_lid_cavity_vtk"
  static std::string dir_from_path(const std::string &path) {
    const std::filesystem::path p(path);
    return (p.parent_path() / (p.stem().string() + "_vtk")).string();
  }

  void run() {
    while (true) {
      std::vector<char> chunk;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [this] { return stop_ || !queue_.empty(); });

        if (queue_.empty() && stop_) {
          break;
        }

        chunk = std::move(queue_.front());
        queue_.pop();
      }

      if (!have_dims_) {
        parse_header(chunk);
      } else {
        write_frame(chunk);
      }
    }
  }

  // Primo chunk atteso: 2 o 3 int32 (nx,ny[,nz]).
  void parse_header(const std::vector<char> &data) {
    const std::size_t n_ints = data.size() / sizeof(int32_t);
    if (n_ints != 2 && n_ints != 3) {
      std::cerr << "[VtkWriter] header inatteso (" << data.size()
                << " byte), ignorato" << std::endl;
      return;
    }

    int32_t nx32 = 0, ny32 = 0, nz32 = 1;
    std::memcpy(&nx32, data.data(), sizeof(int32_t));
    std::memcpy(&ny32, data.data() + sizeof(int32_t), sizeof(int32_t));
    if (n_ints == 3) {
      std::memcpy(&nz32, data.data() + 2 * sizeof(int32_t), sizeof(int32_t));
    }

    nx_ = static_cast<std::size_t>(nx32);
    ny_ = static_cast<std::size_t>(ny32);
    nz_ = static_cast<std::size_t>(nz32);
    have_dims_ = true;
  }

  // Frame successivi: array di float (velocity magnitude), uno per nodo,
  // nello stesso ordine usato da Grid::scalar_index. Ogni frame diventa
  // un file .vti (ImageData XML) con un unico campo scalare in PointData.
  void write_frame(const std::vector<char> &data) {
    const std::size_t n_points = nx_ * ny_ * nz_;
    const std::size_t expected_bytes = n_points * sizeof(float);
    if (data.size() != expected_bytes) {
      std::cerr << "[VtkWriter] frame " << frame_count_ << ": ricevuti "
                << data.size() << " byte, attesi " << expected_bytes
                << ", frame ignorato" << std::endl;
      return;
    }

    std::vector<float> values(n_points);
    std::memcpy(values.data(), data.data(), data.size());

    std::ostringstream fname;
    fname << basename_ << "_" << std::setw(5) << std::setfill('0')
          << frame_count_ << ".vti";
    const std::filesystem::path frame_path =
        std::filesystem::path(output_dir_) / fname.str();

    std::ofstream out(frame_path);
    if (!out.is_open()) {
      std::cerr << "[VtkWriter] impossibile aprire " << frame_path << std::endl;
      return;
    }

    const std::string extent = "0 " + std::to_string(nx_ - 1) + " 0 " +
                               std::to_string(ny_ - 1) + " 0 " +
                               std::to_string(nz_ - 1);

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"ImageData\" version=\"1.0\" "
           "byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "  <ImageData WholeExtent=\"" << extent
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n"
        << "    <Piece Extent=\"" << extent << "\">\n"
        << "      <PointData Scalars=\"velocity_magnitude\">\n"
        << "        <DataArray type=\"Float32\" "
           "Name=\"velocity_magnitude\" format=\"ascii\">\n";

    for (std::size_t i = 0; i < n_points; ++i) {
      out << values[i] << " ";
    }

    out << "\n        </DataArray>\n"
        << "      </PointData>\n"
        << "    </Piece>\n"
        << "  </ImageData>\n"
        << "</VTKFile>\n";
    out.close();

    frame_files_.push_back(fname.str());
    write_pvd();
    ++frame_count_;
  }

  // Riscrive il .pvd (piccolo XML) con tutti i frame prodotti finora, cosi'
  // la serie e' apribile in ParaView anche a simulazione ancora in corso.
  void write_pvd() const {
    const std::filesystem::path pvd_path =
        std::filesystem::path(output_dir_) / (basename_ + ".pvd");

    std::ofstream out(pvd_path);
    if (!out.is_open())
      return;

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n"
        << "  <Collection>\n";

    for (std::size_t i = 0; i < frame_files_.size(); ++i) {
      out << "    <DataSet timestep=\"" << i << "\" group=\"\" part=\"0\" "
          << "file=\"" << frame_files_[i] << "\"/>\n";
    }

    out << "  </Collection>\n"
        << "</VTKFile>\n";
  }

  const std::string output_dir_, basename_;

  bool have_dims_;
  std::size_t nx_, ny_, nz_;
  std::size_t frame_count_;
  std::vector<std::string> frame_files_;

  std::thread worker_;
  std::queue<std::vector<char>> queue_;
  std::mutex mtx_;
  std::condition_variable cv_;
  bool stop_;
};

} // namespace lbm

#endif // __LBM_SIM_DATA_VTK_WRITER
