#include "QSF/qsf.h"
#include "cxxopts.h"
#include <filesystem>
namespace fs = std::filesystem;

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

int main(int argc, char* argv[])
{
	cxxopts::Options options("nitrogen-3e", "3e simulations of nitrogen");
	options.add_options()
		("g,gaussian", "Start from gaussian") // a bool parameter
		("n,nodes", "Number of nodes", cxxopts::value<ind>()->default_value("1024"))
		("f,field", "Field value", cxxopts::value<double>()->default_value("0.12"))
		("d,dt", "Set the timedelta", cxxopts::value<double>()->default_value("0.1"))
		("s,soft", "Set the coulomb softener (epsilon)", cxxopts::value<double>()->default_value("0.92"))
		("b,border", "Number of border nodes defined as nCAP=nodes/#", cxxopts::value<ind>()->default_value("4"));

	auto result = options.parse(argc, argv);

#ifdef BATCH
	const ind nodes = 1024;
	constexpr DIMS my_dim = 3_D;
	//The following gives NEcharge = -3.0 and EEcharge=0.5
	const double Ncharge = 3.0 / sqrt(0.5);
	const double Echarge = -sqrt(0.5);
	const double NEsoft = 1.02;
	const ind nCAP = nodes / 4;
	const double re_dt = 0.1;
	const double omega = 0.06;
	const double F0 = sqrt(2. / 3.) * 0.12;
	fs::path loc{ std::getenv("SCRATCH") };
	fs::path re_output_dir{ IOUtils::project_name };

	re_output_dir = loc / re_output_dir / fs::path("n_" + std::to_string(nodes)) /
		fs::path("nCAP_" + std::to_string(nCAP)) / fs::path("dt_" + std::to_string(re_dt)) / fs::path("F0_" + std::to_string(F0 / sqrt(2. / 3.)));
#else
#ifdef TEST
	const ind nodes = result["nodes"].as<ind>();
	constexpr DIMS my_dim = 2_D;
	const double Ncharge = result["gaussian"].as<bool>() ? 0.0 : 3.23 / sqrt(0.5);
	const double Echarge = result["gaussian"].as<bool>() ? 0.0 : -sqrt(0.5);
	const double NEsoft = result["soft"].as<double>();
	const ind nCAP = nodes / result["border"].as<ind>();
	const double re_dt = result["dt"].as<double>();
	const double omega = 0.06;
	const double F0 = sqrt(2. / 3.) * result["field"].as<double>();
	fs::path loc{ IOUtils::project_dir };
	fs::path re_output_dir{ IOUtils::results_dir };

	re_output_dir = loc / re_output_dir / fs::path("n_" + std::to_string(nodes)) /
		fs::path("nCAP_" + std::to_string(nCAP)) / fs::path("dt_" + std::to_string(re_dt)) / fs::path("F0_" + std::to_string(F0 / sqrt(2. / 3.)));
#else
	const ind nodes = result["nodes"].as<ind>();
	constexpr DIMS my_dim = 3_D;
	const double Ncharge = 3.0 / sqrt(0.5);
	const double Echarge = -sqrt(0.5);
	const double NEsoft = 1.02;
	const ind nCAP = nodes / result["border"].as<ind>();
	const double re_dt = result["dt"].as<double>();
	const double omega = 0.06;
	const double F0 = sqrt(2. / 3.) * result["field"].as<double>();
	fs::path loc{ IOUtils::project_dir };
	fs::path re_output_dir{ IOUtils::results_dir };

	re_output_dir = loc / re_output_dir;
#endif
#endif 


	fs::path gs_dir = loc / fs::path("groundstates");
	fs::path gs_file = fs::path(std::string("nitrogen_") + std::to_string(nodes) + std::string(".psib0"));
	fs::path re_input_file = gs_dir / gs_file;
	logWarning("target folder: %s", re_output_dir.c_str());
	// return 0;
	QSF::init(re_output_dir, argc, argv);
	if (!MPI::pID)
		std::filesystem::create_directories(gs_dir);

	const double dx = 100.0 / 511.0;
	const double ncycles = 3.0;
	const int log_interval = 20;

	const ind halfcycle_steps = round(pi / omega / re_dt);

	EckhardtSachaInteraction potential{ {
		.Ncharge = Ncharge,
		.Echarge = Echarge,
		.Nsoft = NEsoft,
		.Esoft = NEsoft } };

	if (SHOULD_RUN(MODE::IM)) //if any parameter is passed assume gs
	{
		CAP<CartesianGrid<my_dim>> im_grid{ {dx, nodes}, nCAP };
		auto im_wf = Schrodinger::Spin0{ im_grid, potential };
		auto im_outputs = BufferedBinaryOutputs<
			VALUE<Step, Time>
			, OPERATION<Orthogonalize>
			, OPERATION<AntiSymmetrize>
			, OPERATION<Normalize>
			, AVG<Identity>
			, AVG<PotentialEnergy>
			, AVG<KineticEnergy>
			, ENERGY_TOTAL
			, ENERGY_DIFFERENCE
		>{ {.comp_interval = 1, .log_interval = 20} };

		auto p1 = SplitPropagator<MODE::IM, SplitType, decltype(im_wf)>
		{
			{.dt = 0.3, .max_steps = 1000000, .state_accuracy = 10E-15},
			std::move(im_wf)
		};

		p1.run(im_outputs,
			   [&](const WHEN when, const ind step, const uind pass, auto& wf)
			   {
				   if (when == WHEN::AT_START)
				   {
					   wf.addUsingCoordinateFunction(
						   [](auto... x) -> cxd
						   {
							   return cxd{ (x*...*1.0) * gaussian(0.0, 3.0, x...), 0 };
						   });
					   logUser("wf loaded manually!");
				   }
				   if (when == WHEN::AT_END) wf.save(re_input_file);
			   });

	}




	if (SHOULD_RUN(MODE::RE))
	{
		logUser("About to use %s ", re_input_file.c_str());

		CAP<MultiCartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nCAP };
		using F1 = Field<AXIS::XYZ, ChemPhysEnvelope<ChemPhysPulse>>;
		DipoleCoupling<VelocityGauge, F1> re_coupling
		{
			ChemPhysEnvelope<ChemPhysPulse>{ {
				.F0 = F0,
				.omega = omega,
				.ncycles = ncycles,
				.phase_in_pi_units = 0,
				.delay_in_cycles = 0}}
		};

		auto re_outputs = BufferedBinaryOutputs <
			VALUE<Step, Time>
			, VALUE<F1>
			// , AVG<Identity>
			// , AVG<PotentialEnergy>
			// PROJ<EIGENSTATES, Identity>,
			// , AVG<DERIVATIVE<0, PotentialEnergy>>
			// FLUX<BOX<3>>
			, VALUE<ETA>
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto re_wf = Schrodinger::Spin0{ re_capped_grid, potential, re_coupling };
		auto p2 = SplitPropagator<MODE::RE, SplitType, decltype(re_wf)>{ {.dt = re_dt}, std::move(re_wf) };
		p2.run(re_outputs,
			   [=](const WHEN when, const ind step, const uind pass, auto& wf)
			   {
				   if (when == WHEN::AT_START)
				   {
					   if (MPI::region == 0)
					   {
					   #ifdef BATCH 
						   wf.load(re_input_file);
					   #else
						   logUser("Initing from gaussian");
						   if (result["gaussian"].as<bool>())
							   wf.addUsingCoordinateFunction(
								   [](auto... x) -> cxd
								   {
									   double mom = ((x * 8.0) + ...);
									   return gaussian(2.0, 1.0, x...) * cxd { cos(mom), sin(mom) };
								   });
						   else wf.load(re_input_file);
					   #endif
					   }
				   }

				   if (when == WHEN::DURING)
				   {
				   #ifdef BATCH
					   if (step % halfcycle_steps == 0)
						   wf.save("latest_backup");
				   #else
					   if (step == 1 || step % halfcycle_steps == halfcycle_steps - 1)
					   {
						   static int enter = 0;
						   wf.save("during_" + std::to_string(enter));
						   enter++;
					   }
				   #endif 
				   }
				   if (when == WHEN::AT_END)
				   {
					   wf.save("final_");
					   wf.saveIonizedJoined("final_ionized_joined_P_", { my_dim, REP::P, true, true, true, true, false });
				   }
			   });
	}
	QSF::finalize();
}

