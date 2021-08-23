#include "QSF/qsf.h"

#ifdef SIZE
const ind nodes = SIZE;
constexpr DIMS dim = 2_D;
const double Ncharge = 3.23 / sqrt(0.5);
const double Echarge = -sqrt(0.5);
const double NEsoft = 0.92;
const ind nCAP = nodes / 4;
const double re_dt = 0.1;
const double omega = 0.06;
const double F0 = sqrt(2. / 3.) * 0.15;
std::string re_input{ "im0__aft_xy.psib0" };
#else
const ind nodes = 1024;
constexpr DIMS dim = 3_D;
//The following gives NEcharge = -3.0 and EEcharge=0.5
const double Ncharge = 3.0 / sqrt(0.5);
const double Echarge = -sqrt(0.5);
const double NEsoft = 1.02;
const ind nCAP = nodes / 4;
const double re_dt = 0.1;
const double omega = 0.06;
const double F0 = sqrt(2. / 3.) * 0.12;
std::string re_input{ "im0__aft_xyz.psib0" };
#endif 

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

int main(int argc, char* argv[])
{
	std::string loc{ TEST ? project_dir : scratch_dir };
	std::string dir{ TEST ? results_dir : project_name };
	std::string sep{ "/" };

	auto re_input_file = loc + sep + dir + sep + re_input;
	QSF::init(loc.c_str(), dir.c_str(), argc, argv);

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
		CAP<CartesianGrid<dim>> im_grid{ {dx, nodes}, nCAP };
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
				   if (when == WHEN::AT_END) wf.save("im");
			   });

	}

	if (SHOULD_RUN(MODE::RE))
	{
		logInfo("About to use %s ", re_input_file.c_str());

		CAP<MultiCartesianGrid<dim>> re_capped_grid{ {dx, nodes}, nCAP };
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

				//    #if TEST
				// 	   if (MPI::region == 0)
				// 		   wf.addUsingCoordinateFunction(
				// 			   [](auto... x) -> cxd
				// 			   {
				// 				   double mom = ((x * 8.0) + ...);
				// 				   return gaussian(2.0, 1.0, x...) * cxd { cos(mom), sin(mom) };
				// 			   });
				//    #else 
					   if (MPI::region == 0)  wf.load(re_input_file);
				//    #endif
				   }

				   if (when == WHEN::DURING)
				   {
				   #if TEST
					   if (step == 1 || step % halfcycle_steps == 0)
					   {
						   static int enter = 0;
						   wf.save("during_" + std::to_string(enter));
						   enter++;
					   }
				   #else
					   if (step % halfcycle_steps == 0)
						   wf.save("latest_backup");
				   #endif 
				   }

				   if (when == WHEN::AT_END)
				   {
					   wf.save("final_");
					   wf.saveIonizedJoined("final_ionized_joined_P_", { dim, REP::P, true, true, true, true, false });
				   }
			   });
	}
	QSF::finalize();
}

