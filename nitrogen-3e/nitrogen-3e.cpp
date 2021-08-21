#include "QSF/qsf.h"
#include "cxxopts.h"


cxxopts::Options options(
	"nitrogen-3e",
	"This program calculates momenta distributions of nitrogen in reduced-dim Eckhardt-Sacha model");

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;


int main(int argc, char* argv[])
{
	options.add_options()
		("b,batch", "Whether we run on cluster (changes target directories)",
		 cxxopts::value<bool>()->default_value("false"))
		("r,restore", "Whether to continue from the last backup",
		 cxxopts::value<bool>()->default_value("false"))
		("n,nodes", "Set number of nodes manually",
		 cxxopts::value<ind>())
		("t,dt", "Set number timestep manually",
		 cxxopts::value<double>())
		("f,field", "Set number field strength manually",
		 cxxopts::value<double>())
		("i,input", "File name of the initial state for real-time propagation",
		 cxxopts::value<std::string>()->default_value("im0__aft_xyz.psib0"))
		("h,help", "Print usage");
	auto result = options.parse(argc, argv);
	if (result.count("help"))
	{
		logInfo("%s", options.help().c_str());
		exit(0);
	}

	bool batch = result["batch"].as<bool>();
	std::string loc{ batch ? scratch_dir : project_dir };
	std::string dir{ batch ? project_name : results_dir };
	std::string sep{ "/" };

	QSF::init(loc.c_str(), dir.c_str(), argc, argv);

	auto re_input_file = loc + sep + dir + sep + result["input"].as<std::string>();
	const ind nodes = result.count("nodes") ? result["nodes"].as<ind>() : 1024;
	const ind nCAP = nodes / 4;
	const double dx = 100.0 / 511.0;
	const double ncycles = 3.0;
	const double omega = 0.06;
	const double F0 = sqrt(2. / 3.) * (result.count("field") ? result["field"].as<double>() : 0.12);// / omega;
	const double re_dt = result.count("dt") ? result["dt"].as<double>() : 0.05; //delta step used in real time prop.
	const int log_interval = 20;

	const ind halfcycle_steps = round(pi / omega / re_dt);
	//The following gives NEcharge = -3.0 and EEcharge=0.5
	EckhardtSachaInteraction potential{ {
		.Ncharge = 3.0 / sqrt(0.5),
		.Echarge = -sqrt(0.5),
		.Nsoft = 1.02,
		.Esoft = 1.02 } };

	if (SHOULD_RUN(MODE::IM)) //if any parameter is passed assume gs
	{
		CAP<CartesianGrid<3_D>> im_grid{ {dx, nodes}, nCAP };
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

		CAP<MultiCartesianGrid<3_D>> re_capped_grid{ {dx, nodes}, nCAP };
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
		//    , AVG<Identity>, AVG<PotentialEnergy>, PROJ<EIGENSTATES, Identity>, AVG<DERIVATIVE<0, PotentialEnergy>>, FLUX<BOX<3>>
			, VALUE<ETA>
		>{ {.comp_interval = 1, .log_interval = log_interval} };

		auto re_wf = Schrodinger::Spin0{ re_capped_grid, potential, re_coupling };
		auto p2 = SplitPropagator<MODE::RE, SplitType, decltype(re_wf)>{ {.dt = re_dt}, std::move(re_wf) };
		p2.run(re_outputs,
			   [=](const WHEN when, const ind step, const uind pass, auto& wf)
			   {
				   if (when == WHEN::AT_START)
					   if (MPI::region == 0) wf.load(re_input_file);

				   if (when == WHEN::DURING)
					   if (step % halfcycle_steps == 0)
						   wf.save("latest_backup");

				   if (when == WHEN::AT_END)
				   {
					   wf.save("final_");
					//    wf.saveJoined("final_joined_P_", { 3_D, REP::P, true, true, true, true, false });
					   wf.saveIonizedJoined("final_ionized_joined_P_", { 3_D, REP::P, true, true, true, true, false });
				   }
			   });
	}
	QSF::finalize();
}



/* ::::::::::::::::::::::::::::::: Field types :::::::::::::::::::::::::::::::::
	* Each Field may take many pulse types, each of which can be envoloped or not.
	* The first template parameter determines the axis on which the field will act.
	* Fields are automatically converted to vector potential in case of gauges that
	* require vector potentials in their definition, eg. VelocityGauge */

/* :::::::::::::: Generic Hamiltonian type parametrized by fields  :::::::::::::
* The Hamiltonian is parametrized by Kinetic energy, static Potential
* and Coupling - in the gauge of choice. As the Hamiltonian determines the system,
* the wavefunction is generated based on the HamiltonianOperator configuration.*/
/* :::::::::::::::::::::::::: Wavefunction ::::::::::::::::::::::::::::::::::::: */

/* :::::::::::::::::::::: Split operator type & base  ::::::::::::::::::::::::::
* The evolution is performed using the generalizedsplit-operator technique.
* The type can be MultiproductSum or ... both of which require a base to
* be defined. Base can be either kinetic-part or potential-part first TVT/VTV. */
/* ::::::::::::::::::::::::::::: SplitPropagator ::::::::::::::::::::::::::::::::::::
* The Hamiltonian, together with the Split operator type and base determine the
* propagator used in evolution. Moreover we can have evolution in imaginary or real
* time determined by MODE: IM for imaginary time and RE for real time propagation. */

/* ::::::::::::::::::::: Quantities/Computations/Outputs :::::::::::::::::::::::
* Both IM and RE modes need to do some computations, though vastly different ones.
* In RE mode we compute many quantities and we can postpone the moment of MPI reduce
* procedures, while in IM mode we have few quantities which we cannot postpone.
* One could also argue that the quantities of interest in IM mode could be fixed.
* For this reason we predefined them under the name STD_outputsIM. Notice, that the
* computations may be dependent on the coupling, base Hamiltonian or state, thus we
* need to define them as templates, parametrized by T, V, COUPLING and PASS which
* loosly speaking refers to the index of the state currently computed.
* We can also choose quantities which will not be saved in the file (AUXILLARY_VALUES)*/
/*:::::::::::::::::::::::::: Dumps(wf, potential) ::::::::::::::::::::::::::::::
*Sometimes we also need to dump the values of a tensor such as the wavefunction
* or potential at some moments of the computation.We can perform them BEFORE, AFTER
AT_STEP or AT_INTERVAL */



// template <class SingleOrMulti>
// using BaseGrid = CartesianGrid<DIMS::D3, CAP, SingleOrMulti, MPI::Slices>;

// template <class SingleOrMulti, class ... Fields>
// using wf = Schrodinger::Spin0<BaseGrid<SingleOrMulti>, CoulombInteraction, DipoleCoupling<VelocityGauge, Fields...>>;

// using RE_Prop = SplitPropagator<MODE::RE, SplitType, wf<MPI::Multi, F1>>;




// using F1 = Field<AXIS::X, SinPulse>;
// F1 f1
// {
// 	SinPulse{{.F0 = 0.1, .omega = 0.06, .ncycles = 4.0, .FWHM_percent = 0.5, .phase_in_pi_units = 0, .delay_in_cycles = 0}}
// };

// using F2 = Field <AXIS::Y, GaussianEnvelope<SinPulse>>;
// F2 f2
// {
// 	GaussianEnvelope<SinPulse>{ {.F0 = 0.1, .omega = 0.06, .ncycles = 4.0, .FWHM_percent = 0.5, .phase_in_pi_units = 0, .delay_in_cycles = 0}}
// };

// DipoleCoupling<VelocityGauge, F1, F2> re_coupling
// { f1,f2 };


// using wf = Schrodinger::Spin0{ grid_2, potential_2, coupling_2 };

// // <grid_2, CoulombInteraction, DipoleCoupling<VelocityGauge, F1, F2>>;
// using outputs_2 = BufferedBinaryOutputs<
// 	VALUE<Step, Time>
// 	, VALUE<F1, F2>
// >;
// SplitPropagator<MODE::RE, SplitType, wf<MPI::Multi, F1>>

	// auto p2 = RE_Prop{ wf<MPI::Multi,F1>{
	// 	BaseGrid(32,0.5),
	// 	{
	// 		.Ncharge = 1,
	// 		.Echarge = 1,
	// 		.Nsoft = 1,
	// 		.Esoft = 1
	// 	},
	// 	{

	// 	}
	// 	} };
// auto r1 = Repeat<up_to_t<0>, IM_Prop>(
// 	[=](ind step, uind PASS, const auto& wf)
// 	{
// 		if (step % 100 == 0)
// 			logInfo("I'm working at step %td!", ind(step));
// 	}
// );
// r1.run();


		// typename WF::Coupling,/* Field/Potential values and energies */
		// CO_AVG <typename WF::H::T>, CO_AVG <typename WF::H::V>,
	// AVG<AXIS::XYZ, KinEnergy>,
	// AVG<AXIS::XYZ, PotEnergy>>;
	// AVG<AXIS::X, Derivative<PotEnergy>,
	// PROJ<EIGENSTATES, Identity>
	// CO_AVG <IdentityOperator<>>, /* Norm */
	// CO_DIR_AVG <PotentialDerivativeX<typename WF::H::V>>, /* Derivative of potential */	
	// CO_PROJ<IdentityOperator<>, PASSES>, /* Projections on eigenstates */
	// STD_fluxes<typename WF::Coupling>>; /* Fluxes */
	// AUXILLARY_VALUES<ETAOperator>>; /* Estimated end of computation */


	// template <ind PASS>
	// using dumpsIM = Dumps < AFTER<DUMP_PSI<REP::X, DIMS::ALL>>>;
	// template <ind PASS>
	// using dumpsRE = Dumps <
	// 	AFTER<DUMP_PSI<REP::X, DIMS::ALL>>,
	// 	AT_INTERVAL <50, DUMP_PSI<REP::X, DIMS::ALL>>,
	// 	AT_STEP <1, DUMP_PSI<REP::X, DIMS::ALL>>>;

	/* :::::::::::::::::::::::::: Optimizations ::::::::::::::::::::::::::::::::::::::::: */


/* ::::::::::::::::::::::::::: Main configuration :::::::::::::::::::::::::::::::::::
* So far, the alias names used in all "using ... =" statements have been arbitrary.
* They were defined as template aliases on purpose, in order to customize them now.
* PROJECT_CONFIG merges all the above configuration and tells the program how to run. */
// using PROJECT_CONFIG = Run<
	// Routine<up_to_t<0>, IM_Prop, STD_outputsIM, dumpsIM, opt>,
	// Routine<up_to_t<0>, RE_Prop, outputsRE, dumpsRE, opt>>;

// template <ind PASS>
// using conf = CONFIG<INIT, COMP, DUMP>

