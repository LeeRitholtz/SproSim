#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "sprosim/interfaces/IParticle.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/CoffeeBed.h"
//#include "sprosim/WaterFlow.h"

namespace sprosim {
namespace visualization {

/**
 * @brief VTK exporter for ParaView visualization
 *
 * Exports simulation data to VTK format files that can be loaded into ParaView
 * for scientific visualization and analysis. Supports both single timestep
 * exports and time series collections.
 */
class VTKExporter {
public:
    VTKExporter() = default;
    ~VTKExporter() = default;

    /**
     * @brief Export coffee particles to VTK unstructured grid format
     * @param particles Vector of coffee particles to export
     * @param filename Output filename (should end in .vtu)
     * @param timestep Optional timestep number for time series
     */
    void export_particles(const std::vector<std::shared_ptr<ICoffeeParticle>>& particles,
                         const std::string& filename,
                         int timestep = -1);

    /**
     * @brief Export flow field to VTK structured grid format
     * @param flow Water flow field to export
     * @param filename Output filename (should end in .vts)
     * @param timestep Optional timestep number for time series
     */
    void export_flow_field(const std::shared_ptr<IWaterFlow>& flow,
                          const std::string& filename,
                          int timestep = -1);

    /**
     * @brief Export complete simulation state (particles + flow)
     * @param bed Coffee bed containing particles
     * @param flow Water flow field
     * @param base_filename Base name for output files
     * @param timestep Timestep number
     */
    void export_simulation_state(const std::shared_ptr<CoffeeBed>& bed,
                                const std::shared_ptr<IWaterFlow>& flow,
                                const std::string& base_filename,
                                int timestep);

    /**
     * @brief Start a time series export sequence
     * @param collection_filename Name for the collection file (.pvd)
     * @param base_filename Base name for individual timestep files
     */
    void start_time_series(const std::string& collection_filename,
                          const std::string& base_filename);

    /**
     * @brief Add timestep to current time series
     * @param bed Coffee bed state
     * @param flow Flow field state
     * @param time Physical time value
     * @param timestep Timestep number
     */
    void add_timestep_to_series(const std::shared_ptr<CoffeeBed>& bed,
                               const std::shared_ptr<IWaterFlow>& flow,
                               double time,
                               int timestep);

    /**
     * @brief Finalize and close time series collection
     */
    void finalize_time_series();

    /**
     * @brief Set output directory for all exports
     * @param directory Path to output directory (created if doesn't exist)
     */
    void set_output_directory(const std::string& directory);

private:
    struct TimeSeriesEntry {
        std::string filename;
        double time;
        int timestep;
    };

    // Helper methods for VTK file writing
    void write_vtu_header(std::ofstream& file, size_t num_points, size_t num_cells);
    void write_vtu_points(std::ofstream& file,
                         const std::vector<std::shared_ptr<ICoffeeParticle>>& particles);
    void write_vtu_cells(std::ofstream& file, size_t num_particles);
    void write_vtu_point_data(std::ofstream& file,
                             const std::vector<std::shared_ptr<ICoffeeParticle>>& particles);
    void write_vtu_footer(std::ofstream& file);

    void write_vts_header(std::ofstream& file, size_t nx, size_t ny, size_t nz,
                         double dx, double dy, double dz);
    void write_vts_points(std::ofstream& file, size_t nx, size_t ny, size_t nz,
                         double dx, double dy, double dz);
    void write_vts_point_data(std::ofstream& file,
                             const std::shared_ptr<IWaterFlow>& flow);
    void write_vts_footer(std::ofstream& file);

    void write_pvd_header(std::ofstream& file);
    void write_pvd_timesteps(std::ofstream& file);
    void write_pvd_footer(std::ofstream& file);

    // Utility methods
    std::string get_full_path(const std::string& filename) const;
    void ensure_output_directory_exists() const;
    std::string encode_base64(const std::vector<float>& data) const;
    std::string encode_base64(const std::vector<int>& data) const;

    // Member variables
    std::string output_directory_ = "./vtk_output";
    std::string time_series_collection_file_;
    std::string time_series_base_filename_;
    std::vector<TimeSeriesEntry> time_series_entries_;
    bool time_series_active_ = false;

    // VTK format constants
    static constexpr const char* VTK_XML_VERSION = "1.0";
    static constexpr const char* VTK_BYTE_ORDER = "LittleEndian";
    static constexpr const char* VTK_HEADER_TYPE = "UInt64";
};

} // namespace visualization
} // namespace sprosim
