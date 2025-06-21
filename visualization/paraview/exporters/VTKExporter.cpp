#include "VTKExporter.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace sprosim {
namespace visualization {

void VTKExporter::export_particles(const std::vector<std::shared_ptr<ICoffeeParticle>>& particles,
                                   const std::string& filename, int timestep) {
    ensure_output_directory_exists();
    std::string full_path = get_full_path(filename);

    std::ofstream file(full_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + full_path);
    }

    // Write VTU header
    write_vtu_header(file, particles.size(), particles.size());

    // Write particle positions as points
    write_vtu_points(file, particles);

    // Write cells (each particle is a vertex cell)
    write_vtu_cells(file, particles.size());

    // Write particle data (extraction state, concentration, size)
    write_vtu_point_data(file, particles);

    // Write footer
    write_vtu_footer(file);

    file.close();

    std::cout << "Exported " << particles.size() << " particles to " << filename << std::endl;
}

void VTKExporter::export_flow_field(const std::shared_ptr<IWaterFlow>& flow,
                                    const std::string& filename, int timestep) {
    auto [nx, ny] = flow->get_grid_dimensions();
    double cell_size = flow->get_cell_size();

    ensure_output_directory_exists();
    std::string full_path = get_full_path(filename);

    std::ofstream file(full_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + full_path);
    }

    // Write VTS header for structured grid
    write_vts_header(file, nx, ny, 1, cell_size, cell_size, cell_size);

    // Write grid points
    write_vts_points(file, nx, ny, 1, cell_size, cell_size, cell_size);

    // Write flow field data
    write_vts_point_data(file, flow);

    // Write footer
    write_vts_footer(file);

    file.close();

    std::cout << "Exported " << nx << "x" << ny << " flow field to " << filename << std::endl;
}

void VTKExporter::export_simulation_state(const std::shared_ptr<CoffeeBed>& bed,
                                          const std::shared_ptr<IWaterFlow>& flow,
                                          const std::string& base_filename, int timestep) {
    std::stringstream ss;
    ss << base_filename << "_particles_" << std::setfill('0') << std::setw(6) << timestep << ".vtu";
    export_particles(bed->get_particles(), ss.str(), timestep);

    ss.str("");
    ss << base_filename << "_flow_" << std::setfill('0') << std::setw(6) << timestep << ".vts";
    export_flow_field(flow, ss.str(), timestep);
}

void VTKExporter::start_time_series(const std::string& collection_filename,
                                    const std::string& base_filename) {
    time_series_collection_file_ = collection_filename;
    time_series_base_filename_ = base_filename;
    time_series_entries_.clear();
    time_series_active_ = true;

    std::cout << "Started time series export: " << collection_filename << std::endl;
}

void VTKExporter::add_timestep_to_series(const std::shared_ptr<CoffeeBed>& bed,
                                         const std::shared_ptr<IWaterFlow>& flow, double time,
                                         int timestep) {
    if (!time_series_active_) {
        throw std::runtime_error("Time series not started. Call start_time_series() first.");
    }

    // Export particles
    std::stringstream particles_filename;
    particles_filename << time_series_base_filename_ << "_particles_" << std::setfill('0')
                       << std::setw(6) << timestep << ".vtu";
    export_particles(bed->get_particles(), particles_filename.str(), timestep);

    // Export flow field
    std::stringstream flow_filename;
    flow_filename << time_series_base_filename_ << "_flow_" << std::setfill('0') << std::setw(6)
                  << timestep << ".vts";
    export_flow_field(flow, flow_filename.str(), timestep);

    // Add to time series entries
    TimeSeriesEntry entry;
    entry.filename = particles_filename.str();
    entry.time = time;
    entry.timestep = timestep;
    time_series_entries_.push_back(entry);
}

void VTKExporter::finalize_time_series() {
    if (!time_series_active_) {
        return;
    }

    ensure_output_directory_exists();
    std::string full_path = get_full_path(time_series_collection_file_);

    std::ofstream file(full_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open collection file: " + full_path);
    }

    write_pvd_header(file);
    write_pvd_timesteps(file);
    write_pvd_footer(file);

    file.close();
    time_series_active_ = false;

    std::cout << "Finalized time series with " << time_series_entries_.size() << " timesteps in "
              << time_series_collection_file_ << std::endl;
}

void VTKExporter::set_output_directory(const std::string& directory) {
    output_directory_ = directory;
}

// Private helper methods

void VTKExporter::write_vtu_header(std::ofstream& file, size_t num_points, size_t num_cells) {
    file << "<?xml version=\"" << VTK_XML_VERSION << "\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"" << VTK_BYTE_ORDER
         << "\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells
         << "\">\n";
}

void VTKExporter::write_vtu_points(std::ofstream& file,
                                   const std::vector<std::shared_ptr<ICoffeeParticle>>& particles) {
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& particle : particles) {
        auto [x, y] = particle->get_position();
        // Check if particle supports 3D positioning
        double z = 0.0;
        // For now, use y coordinate as depth proxy for realistic 3D visualization
        z = y / 2.0; // Simple mapping: y position becomes depth
        file << "          " << std::fixed << std::setprecision(6) << x << " " << y << " " << z
             << "\n";
    }

    file << "        </DataArray>\n";
    file << "      </Points>\n";
}

void VTKExporter::write_vtu_cells(std::ofstream& file, size_t num_particles) {
    file << "      <Cells>\n";

    // Connectivity (each particle is a single vertex)
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_particles; i++) {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= num_particles; i++) {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";

    // Cell types (1 = vertex)
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_particles; i++) {
        file << "          1\n";
    }
    file << "        </DataArray>\n";

    file << "      </Cells>\n";
}

void VTKExporter::write_vtu_point_data(
    std::ofstream& file, const std::vector<std::shared_ptr<ICoffeeParticle>>& particles) {
    file << "      <PointData>\n";

    // Extraction state
    file << "        <DataArray type=\"Float32\" Name=\"extraction_state\" format=\"ascii\">\n";
    for (const auto& particle : particles) {
        file << "          " << std::fixed << std::setprecision(6)
             << particle->get_extraction_state() << "\n";
    }
    file << "        </DataArray>\n";

    // Concentration
    file << "        <DataArray type=\"Float32\" Name=\"concentration\" format=\"ascii\">\n";
    for (const auto& particle : particles) {
        file << "          " << std::fixed << std::setprecision(6) << particle->get_concentration()
             << "\n";
    }
    file << "        </DataArray>\n";

    // Particle size
    file << "        <DataArray type=\"Float32\" Name=\"particle_size\" format=\"ascii\">\n";
    for (const auto& particle : particles) {
        file << "          " << std::fixed << std::setprecision(6) << particle->get_size() << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
}

void VTKExporter::write_vtu_footer(std::ofstream& file) {
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}

void VTKExporter::write_vts_header(std::ofstream& file, size_t nx, size_t ny, size_t nz, double dx,
                                   double dy, double dz) {
    file << "<?xml version=\"" << VTK_XML_VERSION << "\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"" << VTK_BYTE_ORDER
         << "\">\n";
    file << "  <StructuredGrid WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 "
         << (nz - 1) << "\">\n";
    file << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1)
         << "\">\n";
}

void VTKExporter::write_vts_points(std::ofstream& file, size_t nx, size_t ny, size_t nz, double dx,
                                   double dy, double dz) {
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;
                file << "          " << std::fixed << std::setprecision(6) << x << " " << y << " "
                     << z << "\n";
            }
        }
    }

    file << "        </DataArray>\n";
    file << "      </Points>\n";
}

void VTKExporter::write_vts_point_data(std::ofstream& file,
                                       const std::shared_ptr<IWaterFlow>& flow) {
    auto [nx, ny] = flow->get_grid_dimensions();

    file << "      <PointData>\n";

    // Velocity field
    file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            auto [vx, vy] = flow->get_velocity(i, j);
            file << "          " << std::fixed << std::setprecision(6) << vx << " " << vy
                 << " 0.0\n";
        }
    }
    file << "        </DataArray>\n";

    // Velocity magnitude
    file << "        <DataArray type=\"Float32\" Name=\"velocity_magnitude\" format=\"ascii\">\n";
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            auto [vx, vy] = flow->get_velocity(i, j);
            double magnitude = std::sqrt(vx * vx + vy * vy);
            file << "          " << std::fixed << std::setprecision(6) << magnitude << "\n";
        }
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
}

void VTKExporter::write_vts_footer(std::ofstream& file) {
    file << "    </Piece>\n";
    file << "  </StructuredGrid>\n";
    file << "</VTKFile>\n";
}

void VTKExporter::write_pvd_header(std::ofstream& file) {
    file << "<?xml version=\"" << VTK_XML_VERSION << "\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"" << VTK_BYTE_ORDER
         << "\">\n";
    file << "  <Collection>\n";
}

void VTKExporter::write_pvd_timesteps(std::ofstream& file) {
    for (const auto& entry : time_series_entries_) {
        file << "    <DataSet timestep=\"" << std::fixed << std::setprecision(6) << entry.time
             << "\" group=\"\" part=\"0\" file=\"" << entry.filename << "\"/>\n";
    }
}

void VTKExporter::write_pvd_footer(std::ofstream& file) {
    file << "  </Collection>\n";
    file << "</VTKFile>\n";
}

std::string VTKExporter::get_full_path(const std::string& filename) const {
    return output_directory_ + "/" + filename;
}

void VTKExporter::ensure_output_directory_exists() const {
    try {
        std::filesystem::create_directories(output_directory_);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to create output directory: " + output_directory_ +
                                 ". Error: " + e.what());
    }
}

std::string VTKExporter::encode_base64(const std::vector<float>& data) const {
    // Simple ASCII output for now - base64 encoding can be added later for efficiency
    return "";
}

std::string VTKExporter::encode_base64(const std::vector<int>& data) const {
    // Simple ASCII output for now - base64 encoding can be added later for efficiency
    return "";
}

} // namespace visualization
} // namespace sprosim