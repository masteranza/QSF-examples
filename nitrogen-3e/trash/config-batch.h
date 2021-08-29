constexpr ind nodes = 1024;
constexpr DIMS my_dim = 3_D;
//The following gives NEcharge = -3.0 and EEcharge=0.5
constexpr double Ncharge = 3.0 / sqrt(0.5);
constexpr double Echarge = -sqrt(0.5);
constexpr double NEsoft = 1.02;
constexpr ind nCAP = nodes / 4;
constexpr double re_dt = 0.1;
constexpr double omega = 0.06;
constexpr double delay_in_cycles = 0.0;
const double phase_in_pi_units = result["phase"].as<double>();
const double F0 = sqrt(2. / 3.) * result["field"].as<double>();
fs::path loc{ std::getenv("SCRATCH") };
fs::path re_output_dir{ IOUtils::project_name };

re_output_dir = loc / re_output_dir / fs::path("F0_" + std::to_string(result["field"].as<double>())) / fs::path("phase_" + std::to_string(phase_in_pi_units) + "pi");