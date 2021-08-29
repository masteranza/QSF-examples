#include "QSF.h"
#include "cxxopts.hpp"
#include <filesystem>
namespace fs = std::filesystem;

cxxopts::Options options("nitrogen-3e", "3e simulations of nitrogen");

cxxopts::ParseResult getOpts(const int argc, char* argv[])
{
	options.add_options()
		("help", "Print help")
		("n,nodes", "Number of nodes (default: 1024)",
		 cxxopts::value<ind>()->default_value("1024"))
		("t,dt", "Set the timedelta",
		 cxxopts::value<double>()->default_value("0.1"))
		("s,soft", "Set the coulomb softener (epsilon)",
		 cxxopts::value<double>()->default_value("0.92"))
		("b,border", "Number of border nodes defined as nCAP=nodes/#",
		 cxxopts::value<ind>()->default_value("4"));

	options.add_options("Environment")
		("c,cluster", "Running on cluster (AGH Prometeusz) (default: false)");

	options.add_options("Testing")
		("g,gaussian", "Start from gaussian (default: false)");

	options.add_options("Laser")
		("f,field", "Field value (default: 0.12)",
		 cxxopts::value<double>()->default_value("0.12"))
		("p,phase", "Phase value (default: 0.0)",
		 cxxopts::value<double>()->default_value("0.0"))
		("d,delay", "Delay (default: 0)",
		 cxxopts::value<double>()->default_value("0.0"));

	return options.parse(argc, argv);
}

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

int main(const int argc, char* argv[])
{
	auto result = getOpts(argc, argv);
	const bool cluster = result["cluster"].as<bool>();
	const bool testing = result["gaussian"].as<bool>();
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

	fs::path loc{ cluster ? std::getenv("SCRATCH") : IOUtils::project_dir };
	fs::path re_output_dir{ cluster ? IOUtils::project_name : IOUtils::results_dir };
	re_output_dir = loc / re_output_dir / (testing ? "test" :
										   fs::path("F0_" + std::to_string(result["field"].as<double>()))
										   / fs::path("phase_" + std::to_string(phase_in_pi_units) + "pi"));


	fs::path gs_dir = loc / fs::path("groundstates");
	fs::path gs_file = fs::path(std::string("nitrogen_") + std::to_string(nodes) + std::string(".psib0"));
	fs::path re_input_file = gs_dir / gs_file;

	QSF::init(re_output_dir, argc, argv);
	logImportant("cluster: %d testing: %d", cluster, testing);
	if (!MPI::pID)
	{
		if (result.count("help"))
		{
			std::cout << options.help({ "", "Testing", "Laser" }) << std::endl;
			QSF::finalize();
			exit(0);
		}
		std::filesystem::create_directories(gs_dir);
	}

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
				.phase_in_pi_units = phase_in_pi_units,
				.delay_in_cycles = delay_in_cycles}}
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
		if (testing)
			p2.run(re_outputs, [=](const WHEN when, const ind step, const uind pass, auto& wf)
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0))
					   {
						   wf.addUsingCoordinateFunction(
							   [](auto... x) -> cxd
							   {
								   double mom = ((x * 8.0) + ...);
								   return gaussian(2.0, 1.0, x...) * cxd { cos(mom), sin(mom) };
							   });
					   }
					   if (when == WHEN::DURING)
					   {

						   if (step == 1 || step % halfcycle_steps == halfcycle_steps - 1)
						   {
							   static int enter = 0;
							   wf.save("during_" + std::to_string(enter));
							   enter++;
						   }
					   }
					   if (when == WHEN::AT_END)
					   {
						   wf.save("final_");
						   wf.saveIonizedJoined("final_ionized_joined_P_", { my_dim, REP::P, true, true, true, true, false });
					   }
				   });
		else
			p2.run(re_outputs,
				   [=](const WHEN when, const ind step, const uind pass, auto& wf)
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0))
						   wf.load(re_input_file);

					   if (when == WHEN::DURING)
					   {
						   if (step % halfcycle_steps == 0)
							   wf.save("latest_backup");
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

