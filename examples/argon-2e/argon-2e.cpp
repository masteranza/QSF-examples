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
		("post-core-cut", "Set the range from core to cut after propagation [a.u.]",
		 cxxopts::value<double>()->default_value("50.0"));

   // ("t,dt", "Set the timedelta [a.u.]",
   //  cxxopts::value<double>()->default_value("0.2"))
   // ("s,soft", "Set the coulomb softener (epsilon) [length^2]=[a.u.^2]",
   //  cxxopts::value<double>()->default_value("0.92"))
   // ("b,border", "Number of border nodes defined as nCAP=nodes/# [positive integer]",
   //  cxxopts::value<ind>()->default_value("4"));

	options.add_options("Environment")
		("r,remote", "Running on remote cluster (AGH Prometeusz) (default: false)")
		("k,continue", "Continue calculations from the latest backup (default: false)");

	options.add_options("Testing")
		("g,gaussian", "Start from gaussian wavepacket (default: false)")
		("m,momentum", "Initial momentum of gaussian wavepacket",
		 cxxopts::value<double>()->default_value("3.0"));

	options.add_options("Laser")
		("f,field", "Field strength [a.u] value (default: 8*10^13 W/cm2)",
		 cxxopts::value<double>()->default_value("0.0477447629"))
		("p,phase", "Carrier Envelope Phase (CEP) [pi] value",
		 cxxopts::value<double>()->default_value("0.0"))
		("d,delay", "pulse delay [cycles]",
		 cxxopts::value<double>()->default_value("0.0"))
		("o,postdelay", "post-pulse delay [cycles]",
		 cxxopts::value<double>()->default_value("1.0"))
		("w,fwhm", "FWHM [cycles]",
		 cxxopts::value<double>()->default_value("6.0"));

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
	const double testing_momentum = result["momentum"].as<double>();
	logWarning("testing %d momentum %g", testing, testing_momentum);
	const ind nodes = result["nodes"].as<ind>();
	constexpr DIMS my_dim = 2_D;
	const double Ncharge = 2.0;
	const double Echarge = -1.0;
	const double NEsoft = (testing ? 1000000.0 : 2.2);
	const ind nCAP = 32;//nodes / result["border"].as<ind>();
	const double omega = 0.0146978556546; //3100um
	const double FWHM_cycles = result["fwhm"].as<double>();
	const double delay_in_cycles = result["delay"].as<double>();
	const double postdelay_in_cycles = result["postdelay"].as<double>();
	//The value 3.3 is empirical giving a nice smooth gaussian tail tending towards zero 
	const double ncycles = (testing ? 0.1 : round(FWHM_cycles * 3.3));
	const double phase_in_pi_units = result["phase"].as<double>();
	const double field = testing ? 0.0 : result["field"].as<double>();
	const double dx = result["dx"].as<double>();
	const double gsrcut = result["post-core-cut"].as<double>();
	const double re_dt = dx * 0.25;//dx * dx * 0.5;//result["dt"].as<double>();
	const int log_interval = testing ? 20 : 1000;
	// backup interval should be a multiple of log_interval
	const ind ncycle_steps = log_interval * (testing ? 4 : ind(round(twopi / omega / re_dt) / log_interval));

	IO::path output_dir{ remote ? std::getenv("SCRATCH") : IO::project_dir };
	output_dir /= remote ? IO::project_name : IO::results_dir;
	if (testing) output_dir /= "test";

	QSF::init(argc, argv, output_dir);

	if (result.count("help"))
	{
		if (!MPI::pID)
			std::cout << options.help({ "", "Laser", "Testing" }) << std::endl;
		QSF::finalize();
		exit(0);
	}

	EckhardtSachaInteraction potential{ {
		.Ncharge = Ncharge, .Echarge = Echarge,
		.Nsoft = NEsoft, .Esoft = NEsoft } };

	if constexpr MODE_FILTER_OPT(MODE::IM)
	{
		QSF::subdirectory("groundstates");
		CAP<CartesianGrid<my_dim>> im_grid{ {dx, nodes}, nCAP };
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
			   });
	}

	// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// Real-time part :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	if constexpr MODE_FILTER_OPT(MODE::RE)
	{
		QSF::subdirectory(testing ? "testing" : "fwhm_cycles_" + std::to_string(FWHM_cycles));
		// We need to pass absolute path 
		IO::path im_output = IO::root_dir / IO::path("groundstates/" + std::to_string(nodes) + "_" + std::to_string(dx) + "_repX");
	#define MULTIGRID 1
	#ifdef MULTIGRID
		CAP<MultiCartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nodes / 4 };
	#else
		CAP<CartesianGrid<my_dim>> re_capped_grid{ {dx, nodes}, nCAP };
	#endif
		using A1 = VectorPotential<AXIS::XY, GaussianEnvelope<SinPulse>, ConstantPulse>;
		DipoleCoupling<VelocityGauge, A1> re_coupling
		{
			A1{
			GaussianEnvelope<SinPulse>{ {
				.field = field * sqrt(3.) * 0.5 / omega,
				.omega = omega,
				.ncycles = ncycles,
				.FWHM_percent = FWHM_cycles / ncycles,
				.phase_in_pi_units = phase_in_pi_units,
				.delay_in_cycles = delay_in_cycles}},
			ConstantPulse { {
				.field = field,
				.omega = omega,
				.ncycles = postdelay_in_cycles,
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
				   if (testing)
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0))
					   {
						   wf.addUsingCoordinateFunction(
							   [=](auto... x) -> cxd
							   {
								   double mom = ((x * testing_momentum) + ...);
								   return gaussian(0.0, 4.0, x...) * cxd { cos(mom), sin(mom) };
							   });
					   }
					   if (when == WHEN::DURING)
					   {
						   if (step == 1 || step % ncycle_steps == ncycle_steps - 1)
						   {
							   wf.snapshot("_step" + std::to_string(step), DUMP_FORMAT{ .dim = my_dim, .rep = REP::X });
							   wf.snapshot("_step" + std::to_string(step), DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
							   wf.backup(step);
						   }

					   }
					   if (when == WHEN::AT_END)
					   {
						   wf.snapshot("_step" + std::to_string(step), DUMP_FORMAT{ .dim = my_dim, .rep = REP::X });
						   wf.snapshot("_step" + std::to_string(step), DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
						   if (MPI::region == 3) wf.save("momenta", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
						   wf.save("sep_final", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
						   wf.saveIonizedJoined("final", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });

					   }

				   }
				   else
				   {
					   if ((when == WHEN::AT_START) && (MPI::region == 0)) wf.load(im_output);
					   else if (when == WHEN::DURING && (step % ncycle_steps == 0))
					   {
						   wf.backup(step);
						   //    wf.snapshot("X", DUMP_FORMAT{ .dim = my_dim, .rep = REP::X });
						   //    wf.snapshot("P", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
					   }
					   else if (when == WHEN::AT_END)
					   {
						//    wf.save("final");
						//    wf.saveIonizedJoined("final_p", { .dim = my_dim, .rep = REP::P });
						//    wf.orthogonalizeWith(im_output);
					   #ifdef MULTIGRID
						   if (MPI::region == 3) wf.save("momenta", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 1 });
						   if (MPI::region == 3) wf.save("momenta", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 2 });
						   if (MPI::region == 3) wf.save("momenta", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 4 });
						   if (MPI::region == 3) wf.save("momenta", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 8 });
						   wf.saveIonizedJoined("final", DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
					   #else
						   auto name = "n" + std::to_string(nodes) + "_dx" + std::to_string(dx) + "_dt" + std::to_string(re_dt)
							   + "_p" + std::to_string(postdelay_in_cycles) + "_g" + std::to_string(gsrcut);

						   wf.croossOut(gsrcut);
						   wf.save(name, DUMP_FORMAT{ .dim = my_dim, .rep = REP::P });
						   wf.save(name, DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 2 });
						   wf.save(name, DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 4 });
						   wf.save(name, DUMP_FORMAT{ .dim = my_dim, .rep = REP::P, .downscale = 8 });
					   #endif
					   }
				   }

			   }, continu);
	}


	QSF::finalize();
}

