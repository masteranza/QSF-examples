const ind nodes = result["nodes"].as<ind>();
constexpr DIMS my_dim = 2_D;
const double Ncharge = result["gaussian"].as<bool>() ? 0.0 : 3.23 / sqrt(0.5);
const double Echarge = result["gaussian"].as<bool>() ? 0.0 : -sqrt(0.5);
const double NEsoft = result["soft"].as<double>();
const ind nCAP = nodes / result["border"].as<ind>();
const double re_dt = result["dt"].as<double>();
constexpr double omega = 0.06;
constexpr double delay_in_cycles = result["delay"].as<double>();
const double phase_in_pi_units = result["phase"].as<double>();
const double field = sqrt(2. / 3.) * result["field"].as<double>();

fs::path loc{ IOUtils::project_dir };
fs::path re_output_dir{ IOUtils::results_dir };

re_output_dir = loc / re_output_dir / fs::path("n_" + std::to_string(nodes)) /
fs::path("nCAP_" + std::to_string(nCAP)) / fs::path("dt_" + std::to_string(re_dt)) / fs::path("F0_" + std::to_string(field / sqrt(2. / 3.)));