#include "QSF.h"
#include "cxxopts.hpp"
#include <filesystem>

cxxopts::Options options("argon-2e", "2e simulations of nitrogen");

cxxopts::ParseResult getOpts(const int argc, char* argv[])
{
	options.add_options()
		("h,help", "Print help")
		("n,nodes", "Number of nodes [positive integer]",
		 cxxopts::value<ind>()->default_value("1024"))
		("x,dx", "Set the spacedelta [a.u.]",
		 cxxopts::value<double>()->default_value("0.2"))
		("t,dt", "Set the timedelta [a.u.]",
		 cxxopts::value<double>()->default_value("0.2"))
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
		("f,field", "Field strength [a.u] value (default: 8*10^13 W/cm2)",
		 cxxopts::value<double>()->default_value("0.08"))
		("o,omega", "Frequency [au]",
		 cxxopts::value<double>()->default_value("0.094"))
		("p,phase", "Carrier Envelope Phase (CEP) [pi] value",
		 cxxopts::value<double>()->default_value("0.0"))
		("d,delay", "pulse delay [cycles]",
		 cxxopts::value<double>()->default_value("0.0"))
		("w,fwhm", "FWHM [cycles]",
		 cxxopts::value<double>()->default_value("6.0"))
		("c,cycles", "Cycles [number]",
		 cxxopts::value<double>()->default_value("5.0")->implicit_value("5.0"));

	return options.parse(argc, argv);
}

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;
// using TVT = Split3Base<REP::P, REP::X, REP::P>;
// using SplitType = MultiProductSplit<TVT, order>;

int main(const int argc, char* argv[])
{
	using namespace QSF;

	auto result = getOpts(argc, argv);
	const bool remote = result["remote"].as<bool>();
	const bool continu = result["continue"].as<bool>();
	const bool testing = result["gaussian"].as<bool>();
	const bool testing_momentum = result["momentum"].as<double>();
	const ind nodes = result["nodes"].as<ind>();
	constexpr DIMS my_dim = 2_D;
	const double Ncharge = 2.0;
	const double Echarge = -1.0;
	const double NEsoft = 0.6;
	const ind nCAP = 64;//nodes / result["border"].as<ind>();
	const double omega = result["omega"].as<double>();//0.094; //0.023
	const double FWHM_cycles = result["fwhm"].as<double>();
	const double delay_in_cycles = result["delay"].as<double>();
	const double ncycles = 5.0;//result["cycles"].as<double>();
	const double phase_in_pi_units = result["phase"].as<double>();
	const double field = result["field"].as<double>();
	const double dx = result["dx"].as<double>();
	const double re_dt = dx * dx * 0.5;//result["dt"].as<double>();
	const int log_interval = 1000;
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
		.Ncharge = Ncharge, .Echarge = Echarge,
		.Nsoft = NEsoft, .Esoft = NEsoft } };

	if constexpr MODE_FILTER_OPT(MODE::IM)
	{
		QSF::subdirectory("groundstates");
		// CAP<CartesianGrid<my_dim>> im_grid{ {dx, nodes}, nCAP };
		CartesianGrid<my_dim> im_grid{ dx, nodes };
		auto im_wf = Schrodinger::Spin0{ im_grid, potential };
		auto im_outputs = BufferedBinaryOutputs<
			VALUE<Step, Time>
			, OPERATION<Normalize>
			, AVG<Identity>
			, AVG<PotentialEnergy>
			, AVG<KineticEnergy>
			, ENERGY_TOTAL
			, ENERGY_DIFFERENCE
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto p1 = SplitPropagator<MODE::IM, SplitType, decltype(im_wf)>
		{
			{.dt = re_dt, .max_steps = 1000000, .state_accuracy = 10E-15},
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
							   return cxd{ gaussian(0.0, 5.0, x...), 0 };
						   });
					   logUser("wf loaded manually!");
				   }
				   if (when == WHEN::AT_END)
					   wf.save(std::to_string(nodes) + "_" + std::to_string(dx));
			   }, continu);
	}

	// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// Real-time part :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	if constexpr MODE_FILTER_OPT(MODE::RE)
	{
		QSF::subdirectory("n" + std::to_string(nodes) + "_" + std::to_string(dx) + "/F0_" + std::to_string(field) + "/omega_" + std::to_string(omega));
		// We need to pass absolute path 
		IO::path im_output = IO::root_dir / IO::path("groundstates/" + std::to_string(nodes) + "_" + std::to_string(dx));

		// MultiCartesianGrid
		// CAP<MultiCartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nCAP };
		CAP<CartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nCAP };
		// using F1 = Field<AXIS::XY, Sin2Envelope<SinPulse>>;
		using A1 = VectorPotential<AXIS::XY, Sin2Envelope<SinPulse>, ConstantPulse>;
		DipoleCoupling<VelocityGauge, A1> re_coupling
		{
			A1{
			Sin2Envelope<SinPulse>{ {
				.field = field * sqrt(3.) * 0.5 / omega,
				.omega = omega,
				.ncycles = ncycles,
				.phase_in_pi_units = phase_in_pi_units,
				.delay_in_cycles = delay_in_cycles}},
			ConstantPulse { {
				.field = 0.0,
				.omega = omega,
				.ncycles = 3.0,
				.delay_in_cycles = ncycles}}
			}
		};

		auto re_outputs = BufferedBinaryOutputs <
			VALUE<Step, Time>
			, VALUE<A1>
			, AVG<Identity>
			// , AVG<PotentialEnergy>
			// , AVG<KineticEnergy>
			// , AVG<DERIVATIVE<0, PotentialEnergy>>
			// , ZOA_FLUX_2D
			, VALUE<ETA>
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto re_wf = Schrodinger::Spin0{ re_capped_grid, potential, re_coupling };
		auto p2 = SplitPropagator<MODE::RE, SplitType, decltype(re_wf)>{ {.dt = re_dt}, std::move(re_wf) };

		p2.run(re_outputs,
			   [=](const WHEN when, const ind step, const uind pass, auto& wf)
			   {
				   if ((when == WHEN::AT_START) && (MPI::region == 0))
					   wf.load(im_output);
				   else if (when == WHEN::AT_END)
				   {
					//    wf.save("final");
					//    wf.saveIonizedJoined("final_p", { .dim = my_dim, .rep = REP::P });
					//    wf.orthogonalizeWith(im_output);
					   wf.croossOut(50.0);
					   wf.save("finalX", DUMP_FORMAT{ .dim = my_dim, .rep = REP::X });
					   wf.save("finalP", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
				   }

			   }, continu);
	}

	QSF::finalize();
}

