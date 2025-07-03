#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/CoffeeParticle3D.h"
#include "sprosim/PhysicsSolver.h"
#include "sprosim/WaterFlow.h"

#ifdef SPROSIM_HAS_VTK
#include "paraview/exporters/VTKExporter.h"
#endif

using namespace sprosim;

struct BrewingParameters {
    // Core brewing parameters
    double dose = 18.0;           // grams
    int particles = 350;          // particle count
    double tamping_force = 20.0;  // lbs of pressure
    std::string grind = "medium"; // coarse|medium|fine
    double grind_min = 300.0;     // microns (for custom range)
    double grind_max = 600.0;     // microns (for custom range)
    double pressure = 9.0;        // bar
    double temperature = 95.0;    // °C
    double time = 30.0;           // seconds
    double target_yield = 36.0;   // grams of output
    bool use_3d = false;          // particle type
    bool export_vtk = false;      // VTK export

    // Advanced parameters
    double porosity = 0.4;         // initial bed porosity
    double extraction_rate = 0.05; // rate constant
    int grid_width = 58;           // flow field width
    int grid_height = 30;          // flow field height
    double timestep = 0.01;        // simulation dt
};

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  --dose GRAMS          Coffee dose (default: 18.0)\n"
              << "  --particles N         Number of particles (default: 350)\n"
              << "  --tamping-force LBS   Tamping force in lbs (default: 20.0)\n"
              << "  --grind TYPE          Grind size: coarse|medium|fine (default: medium)\n"
              << "  --grind-range MIN MAX Custom grind range in microns\n"
              << "  --pressure BAR        Brewing pressure (default: 9.0)\n"
              << "  --temperature C       Water temperature (default: 95.0)\n"
              << "  --time SECONDS        Brewing time (default: 30.0)\n"
              << "  --yield GRAMS         Target output yield (default: 36.0)\n"
              << "  --3d                  Use 3D particles\n"
              << "  --vtk                 Export VTK files\n"
              << "  --porosity VALUE      Initial bed porosity (default: 0.4)\n"
              << "  --extraction-rate VAL Extraction rate constant (default: 0.05)\n"
              << "  --grid-size WxH       Flow field grid (default: 58x30)\n"
              << "  --timestep DT         Simulation timestep (default: 0.01)\n"
              << "  --help                Show this help\n";
}

BrewingParameters parse_arguments(int argc, char* argv[]) {
    BrewingParameters params;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        } else if (strcmp(argv[i], "--dose") == 0 && i + 1 < argc) {
            params.dose = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--particles") == 0 && i + 1 < argc) {
            params.particles = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--tamping-force") == 0 && i + 1 < argc) {
            params.tamping_force = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--grind") == 0 && i + 1 < argc) {
            params.grind = argv[++i];
        } else if (strcmp(argv[i], "--grind-range") == 0 && i + 2 < argc) {
            params.grind_min = std::stod(argv[++i]);
            params.grind_max = std::stod(argv[++i]);
            params.grind = "custom";
        } else if (strcmp(argv[i], "--pressure") == 0 && i + 1 < argc) {
            params.pressure = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--temperature") == 0 && i + 1 < argc) {
            params.temperature = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--time") == 0 && i + 1 < argc) {
            params.time = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--yield") == 0 && i + 1 < argc) {
            params.target_yield = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--3d") == 0) {
            params.use_3d = true;
        } else if (strcmp(argv[i], "--vtk") == 0) {
            params.export_vtk = true;
        } else if (strcmp(argv[i], "--porosity") == 0 && i + 1 < argc) {
            params.porosity = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--extraction-rate") == 0 && i + 1 < argc) {
            params.extraction_rate = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--grid-size") == 0 && i + 1 < argc) {
            sscanf(argv[++i], "%dx%d", &params.grid_width, &params.grid_height);
        } else if (strcmp(argv[i], "--timestep") == 0 && i + 1 < argc) {
            params.timestep = std::stod(argv[++i]);
        }
    }

    return params;
}

std::pair<double, double> get_grind_range(const BrewingParameters& params) {
    if (params.grind == "coarse")
        return {600, 1000};
    if (params.grind == "fine")
        return {200, 400};
    if (params.grind == "custom")
        return {params.grind_min, params.grind_max};
    return {300, 600}; // medium
}

class ParameterizableBrewingDemo {
  private:
    BrewingParameters params_;

  public:
    ParameterizableBrewingDemo(const BrewingParameters& params) : params_(params) {}

    void run() {
        print_parameters();

        auto coffee_bed = create_coffee_bed();
        auto water_flow = create_water_flow();
        auto solver = create_physics_solver(coffee_bed, water_flow);

        run_simulation(coffee_bed, water_flow, solver);
        display_results(coffee_bed);
    }

  private:
    void print_parameters() {
        std::cout << "=== SproSim Parameterizable Brewing Demo ===\n";
        std::cout << "Dose: " << params_.dose << "g | ";
        std::cout << "Particles: " << params_.particles << " | ";
        std::cout << "Grind: " << params_.grind << " | ";
        std::cout << "Pressure: " << params_.pressure << " bar\n";
        std::cout << "Temperature: " << params_.temperature << "°C | ";
        std::cout << "Time: " << params_.time << "s | ";
        std::cout << "Yield: " << params_.target_yield << "g | ";
        std::cout << "Tamping: " << params_.tamping_force << " lbs | ";
        std::cout << "Type: " << (params_.use_3d ? "3D" : "2D") << "\n\n";
    }

    std::shared_ptr<CoffeeBed> create_coffee_bed() {
        double portafilter_diameter_mm = 58.0;
        auto bed = std::make_shared<CoffeeBed>(params_.dose, portafilter_diameter_mm);

        // Apply tamping force (convert lbs to Pa for compaction)
        double tamping_pa = params_.tamping_force * 6895.0; // lbs to Pa
        bed->update_compaction(tamping_pa);

        // Create particles with specified grind distribution
        auto [grind_min, grind_max] = get_grind_range(params_);
        double bed_radius = portafilter_diameter_mm / 2000.0; // mm to m

        for (int i = 0; i < params_.particles; i++) {
            // Random position within bed
            double angle = 2.0 * M_PI * i / params_.particles;
            double radius_fraction = sqrt(static_cast<double>(i) / params_.particles);
            double x = bed_radius * radius_fraction * cos(angle);
            double y = bed_radius * radius_fraction * sin(angle);

            // Random particle size within grind range
            double size_microns =
                grind_min + (grind_max - grind_min) * (static_cast<double>(i) / params_.particles);
            double radius = size_microns * 1e-6 / 2.0; // microns to meters, diameter to radius

            if (params_.use_3d) {
                double z = 0.0; // Start at bed bottom
                auto particle = std::make_shared<CoffeeParticle3D>(x, y, z, radius);
                bed->add_particle(particle);
            } else {
                auto particle = std::make_shared<CoffeeParticle>(x, y, radius);
                bed->add_particle(particle);
            }
        }

        return bed;
    }

    std::shared_ptr<WaterFlow> create_water_flow() {
        return std::make_shared<WaterFlow>(params_.grid_width, params_.grid_height);
    }

    std::shared_ptr<PhysicsSolver> create_physics_solver(std::shared_ptr<CoffeeBed> bed,
                                                         std::shared_ptr<WaterFlow> flow) {
        sprosim::Parameters solver_params{.permeability = 1e-12,
                                          .fluid_viscosity = 1e-3,
                                          .extraction_rate = params_.extraction_rate,
                                          .temperature = params_.temperature + 273.15, // °C to K
                                          .inlet_pressure = params_.pressure * 1e5 +
                                                            1e5,  // bar to Pa + atmospheric
                                          .outlet_pressure = 1e5, // atmospheric pressure
                                          .particle_drag = 0.1,
                                          .flow_resistance = 0.2,
                                          .saturation_concentration = 0.15,
                                          .temperature_factor = 0.008};

        return std::make_shared<PhysicsSolver>(bed, flow, solver_params);
    }

    void run_simulation(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<WaterFlow> flow,
                        std::shared_ptr<PhysicsSolver> solver) {

        std::cout << "Running simulation...\n";

#ifdef SPROSIM_HAS_VTK
        std::unique_ptr<sprosim::visualization::VTKExporter> vtk_exporter;
        if (params_.export_vtk) {
            vtk_exporter = std::make_unique<sprosim::visualization::VTKExporter>();
            vtk_exporter->set_output_directory("./brew_output");
            vtk_exporter->start_time_series("brewing.pvd", "timestep");
            std::cout << "VTK export enabled -> ./brew_output/\n";
        }
#endif

        const int total_steps = static_cast<int>(params_.time / params_.timestep);
        const int report_interval = total_steps / 10;

        double current_time = 0.0;
        for (int step = 0; step < total_steps; step++) {
            solver->simulate_step(params_.timestep);
            current_time += params_.timestep;

#ifdef SPROSIM_HAS_VTK
            if (params_.export_vtk && vtk_exporter && step % (total_steps / 50) == 0) {
                vtk_exporter->add_timestep_to_series(bed, flow, current_time, step);
            }
#endif

            if (step % report_interval == 0 || step == total_steps - 1) {
                double progress = 100.0 * step / total_steps;
                double extraction = bed->get_extraction_yield() * 100;
                std::cout << std::fixed << std::setprecision(1);
                std::cout << "Progress: " << std::setw(5) << progress << "% | ";
                std::cout << "Extraction: " << std::setw(5) << extraction << "%\n";
            }
        }

#ifdef SPROSIM_HAS_VTK
        if (params_.export_vtk && vtk_exporter) {
            vtk_exporter->finalize_time_series();
        }
#endif
    }

    void display_results(std::shared_ptr<CoffeeBed> bed) {
        std::cout << "\n=== Brewing Results ===\n";

        double extraction = bed->get_extraction_yield() * 100;
        double tds = bed->get_total_dissolved_solids() * 100;

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Extraction Yield: " << extraction << "%\n";
        std::cout << "Total Dissolved Solids: " << tds << "%\n";
        std::cout << "Final Porosity: " << bed->get_porosity() << "\n";

        // Quality assessment
        std::cout << "\nQuality: ";
        if (extraction < 18.0) {
            std::cout << "UNDER-EXTRACTED (increase time/temperature/fineness)\n";
        } else if (extraction > 24.0) {
            std::cout << "OVER-EXTRACTED (decrease time/temperature/fineness)\n";
        } else {
            std::cout << "WELL-EXTRACTED (balanced)\n";
        }

#ifdef SPROSIM_HAS_VTK
        if (params_.export_vtk) {
            std::cout << "\nVisualization: Open ./brew_output/brewing.pvd in ParaView\n";
        }
#endif
    }
};

int main(int argc, char* argv[]) {
    try {
        auto params = parse_arguments(argc, argv);
        ParameterizableBrewingDemo demo(params);
        demo.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
