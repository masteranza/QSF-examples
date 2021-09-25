#include "QSF.h"
#include "cxxopts.hpp"
#include <filesystem>

cxxopts::Options options("neon-3e", "3e simulations of neon like system");

cxxopts::ParseResult getOpts(const int argc, char* argv[])
{
	options.add_options()
		("h,help", "Print help")
		("n,nodes", "Number of nodes [positive integer]",
		 cxxopts::value<ind>()->default_value("1024"))
		("t,dt", "Set the timedelta [a.u.]",
		 cxxopts::value<double>()->default_value("0.1"))
		("s,soft", "Set the coulomb softener (epsilon) [length^2]=[a.u.^2]",
		 cxxopts::value<double>()->default_value("0.92"))
		("b,border", "Number of border nodes defined as nCAP=nodes/# [positive integer]",
		 cxxopts::value<ind>()->default_value("4"));

	options.add_options("Environment")
		("r,remote", "Running on remote cluster (AGH Prometeusz) (default: false)")
		("k,continue", "Continue calculations from the latest backup (default: false)");

	options.add_options("Testing")
		("g,gaussian", "Start from gaussian wavepacket (default: false)")
		("m,momentum", "Initial momentum of gaussian wavepacket",
		 cxxopts::value<double>()->default_value("3.0"));

	options.add_options("Laser")
		("f,field", "Field strength [a.u] value",
		 cxxopts::value<double>()->default_value("0.12"))
		("p,phase", "Carrier Envelope Phase (CEP) [pi] value",
		 cxxopts::value<double>()->default_value("0.0"))
		("d,delay", "pulse delay [cycles]",
		 cxxopts::value<double>()->default_value("0.0"))
		("c,cycles", "Cycles [number]",
		 cxxopts::value<double>()->default_value("3.0")->implicit_value("3.0"));

	return options.parse(argc, argv);
}

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

int main(const int argc, char* argv[])
{
	using namespace QSF;

	auto result = getOpts(argc, argv);
	const bool remote = result["remote"].as<bool>();
	const bool continu = result["continue"].as<bool>();
	const bool testing = result["gaussian"].as<bool>();
	const bool testing_momentum = result["momentum"].as<double>();
	const ind nodes = result["nodes"].as<ind>();
	constexpr DIMS my_dim = 3_D;
	const ind nCAP = nodes / result["border"].as<ind>();
	const double re_dt = result["dt"].as<double>();
	const double omega = 0.06;
	const double delay_in_cycles = result["delay"].as<double>();
	const double ncycles = result["cycles"].as<double>();
	const double phase_in_pi_units = result["phase"].as<double>();
	const double field = result["field"].as<double>();
	const double dx = 100.0 / 511.0;
	const int log_interval = 20;
	// backup interval should be a multiple of log_interval
	const ind halfcycle_steps = log_interval * ind(round(pi / omega / re_dt) / log_interval);

	IO::path output_dir{ remote ? std::getenv("SCRATCH") : IO::project_dir };
	output_dir /= remote ? IO::project_name : IO::results_dir;
	if (testing) output_dir /= "test";
	QSF::init(argc, argv, output_dir);

	logImportant("remote: %d testing: %d", remote, testing);
	if (result.count("help"))
	{
		if (!MPI::pID)
			std::cout << options.help({ "","Environment", "Testing", "Laser" }) << std::endl;
		QSF::finalize();
		exit(0);
	}

	EckhardtSachaInteraction potential{ {
		.Ncharge = 3.0 / sqrt(1.0), .Echarge = -sqrt(1.0),
		.Nsoft = 0.83, .Esoft = 0.83 } };

	if constexpr MODE_FILTER_OPT(MODE::IM)
	{
		QSF::subdirectory("groundstates");
		CAP<CartesianGrid<my_dim>> im_grid{ {dx, nodes}, nCAP };
		auto im_wf = Schrodinger::Spin0{ im_grid, potential };
		auto im_outputs = BufferedBinaryOutputs<
			VALUE<Step, Time>
			, OPERATION<Orthogonalize>
			, OPERATION<AntiSymmetrize<2_D>>
			, OPERATION<Normalize>
			, AVG<Identity>
			, AVG<PotentialEnergy>
			, AVG<KineticEnergy>
			, ENERGY_TOTAL
			, ENERGY_DIFFERENCE
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto p1 = SplitPropagator<MODE::IM, SplitType, decltype(im_wf)>
		{
			PropagatorBase{.dt = 0.3, .max_steps = 1000000, .state_accuracy = 10E-15},
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
							   return cxd{ gaussian(0.0, 3.0, x...), 0 };
						   });
					   logUser("wf loaded manually!");
				   }
				   if (when == WHEN::AT_END)
					   wf.save(std::to_string(nodes));
			   }, continu);
	}

	// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// Real-time part :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	if constexpr MODE_FILTER_OPT(MODE::RE)
	{
		QSF::subdirectory("n" + std::to_string(nodes) + "/F0_" + std::to_string(field) + "/phase_" + std::to_string(phase_in_pi_units) + "pi");
		// We need to pass absolute path 
		IO::path im_output = IO::root_dir / IO::path("groundstates/" + std::to_string(nodes));

		CAP<MultiCartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nCAP };
		using F1 = Field<AXIS::XYZ, ChemPhysEnvelope<ChemPhysPulse>>;
		DipoleCoupling<VelocityGauge, F1> re_coupling
		{
			ChemPhysEnvelope<ChemPhysPulse>{ {
				.field = sqrt(2. / 3.) * field,
				.omega = omega,
				.ncycles = ncycles,
				.phase_in_pi_units = phase_in_pi_units,
				.delay_in_cycles = delay_in_cycles}}
		};

		auto re_outputs = BufferedBinaryOutputs <
			VALUE<Step, Time>
			, VALUE<F1>
			, AVG<Identity>
			// , AVG<PotentialEnergy>
			// , AVG<KineticEnergy>
			// PROJ<EIGENSTATES, Identity>,
			// , AVG<DERIVATIVE<0, PotentialEnergy>>
			// , ZOA_FLUX_3D
			, VALUE<ETA>
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto re_wf = Schrodinger::Spin0{ re_capped_grid, potential, re_coupling };
		auto p2 = SplitPropagator<MODE::RE, SplitType, decltype(re_wf)>{ {.dt = re_dt}, std::move(re_wf) };

		p2.run(re_outputs,
			   [=](const WHEN when, const ind step, const uind pass, auto& wf)
			   {
				   if (!testing)
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0))
						   wf.load(im_output);
					   else if (when == WHEN::DURING && (step % halfcycle_steps == 0))
						   wf.backup(step);
					   else if (when == WHEN::AT_END)
					   {
						   wf.save("final");
						   wf.saveIonizedJoined("final_p", { .dim = my_dim, .rep = REP::P });
					   }
				   }
				   else
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0))
					   {
						   wf.addUsingCoordinateFunction(
							   [=](auto... x) -> cxd
							   {
								   double mom = ((x * testing_momentum) + ...);
								   return gaussian(2.0, 1.0, x...) * cxd { cos(mom), sin(mom) };
							   });
					   }
					   if (when == WHEN::DURING)
					   {
						   if (step == 1 || step % halfcycle_steps == halfcycle_steps - 1)
							   wf.snapshot("_step" + std::to_string(step));

					   }
					   if (when == WHEN::AT_END)
					   {
						   wf.save("final_");
						   wf.saveIonizedJoined("final_p", { .dim = my_dim, .rep = REP::P });
					   }
				   }
			   }, continu);
	}

	QSF::finalize();
}

