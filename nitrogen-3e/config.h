const ind nodes = result["nodes"].as<ind>();
constexpr DIMS my_dim = 3_D;
const double Ncharge = 3.0 / sqrt(0.5);
const double Echarge = -sqrt(0.5);
const double NEsoft = 1.02;
const ind nCAP = nodes / result["border"].as<ind>();
const double re_dt = result["dt"].as<double>();
const double omega = 0.06;
const double delay_in_cycles = result["delay"].as<double>();
const double phase_in_pi_units = result["phase"].as<double>();
const double F0 = sqrt(2. / 3.) * result["field"].as<double>();
fs::path loc{ IOUtils::project_dir };
fs::path re_output_dir{ IOUtils::results_dir };

re_output_dir = loc / re_output_dir;