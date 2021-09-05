#include "QSF.h"

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

int main(int argc, char* argv[])
{
	QSF::init(IOUtils::project_dir / IOUtils::results_dir, argc, argv);
	using im_grid_t = CartesianGrid<2_D>; //2_D == DIMS::D2
	using im_wf_t = Schrodinger::Spin0<im_grid_t, CoulombInteraction>;
	using im_outputs_t = BufferedBinaryOutputs<
		VALUE<Step, Time>
		// , OPERATION<Orthogonalize>
		// , OPERATION<Symmetrize>
		, OPERATION<Normalize>
		, AVG<Identity>
		, AVG<PotentialEnergy>
		, AVG<KineticEnergy>
		, SUM<AVG<PotentialEnergy>, AVG<KineticEnergy>>
		, CHANGE<SUM<AVG<PotentialEnergy>, AVG<KineticEnergy>>>
	>;

	std::string im_output_name = "./Results/im0ELECeDIMd__aft_xyz.psib0";

	{
		auto p1 = SplitPropagator<MODE::IM, SplitType, im_wf_t>{};
		p1.run<im_outputs_t>(
			[&](const WHEN when, const ind step, const uind pass, auto& wf)
			{
				if (when == WHEN::AT_START)
				{
					wf.addUsingCoordinateFunction(
						[](auto... x) -> cxd
						{
							// std::size_t i = 0;
							// double res = 0.0;
							// ((res = x * 1.2, true) || ...); //select first
							// ((i++ == 0 ? (res = x * 0.0, true) : false) || ...); //select i-th
							// return cxd{ cos(res) * gaussian(0.0, 2.0, x...), sin(res) };
							return cxd{ gaussian(0.0, 2.0, x...), 0 };
						});
					logUser("wf loaded manually!");
		// 			// wf.save("ati0_");
				}
		// 		if (when == WHEN::DURING)
		// 		{
		// 			// if (step == 1)wf.save("ati1_");
		// 			// if (step == 50)wf.save("ati50_");
		// 			// if (step == 100)wf.save("ati100_");
		// 			// if (step == 150)wf.save("ati150_");
		// 			// if (step == 200)wf.save("ati200_");
		// 			// if (step == 1000)wf.save("ati1000_");
		// 			// if (step == 10000)wf.save("ati10000_");
		// 		}
		// 		if (when == WHEN::AT_END)
		// 			im_output_name = wf.save("im");

			});
	}
	CAP<CartesianGrid<2_D>> re_capped_grid{ {0.2, 256}, 64 };
	CoulombInteraction re_potential{ {.Ncharge = 1, .Echarge = -1, .Nsoft = 0.0, .Esoft = 0.0 } };

	using F1 = Field<AXIS::XYZ, SinPulse>;
	// using F2 = Field<AXIS::Y, GaussianEnvelope<SinPulse>>;
	DipoleCoupling<VelocityGauge, F1> re_coupling
	{
		SinPulse{{
			// .field = .2, .omega = 0.06,
			.field = .2 * sqrt(2. / 3.), .omega = 0.06,
			.ncycles = 1.0, .FWHM_percent = 0.9,
			.phase_in_pi_units = 0, .delay_in_cycles = 0
			}}
		// ,GaussianEnvelope<SinPulse>{ {
			// .field = 0.1, .omega = 0.06,
			// .ncycles = 4.0, .FWHM_percent = 0.5,
			// .phase_in_pi_units = 0, .delay_in_cycles = 0
			// }}
	};

	auto re_outputs = BufferedBinaryOutputs <
		VALUE<Step, Time>
		, VALUE<F1>
	   // , AVG<Identity>
	   // , AVG<PotentialEnergy>
   // , PROJ<EIGENSTATES, Identity>
	   // , AVG<DERIVATIVE<0, PotentialEnergy>>
	   // , FLUX<BOX<3>>
	   // AUXILLARY_VALUES<ETAOperator>>; /* Estimated end of computation */
	>{ {.comp_interval = 1, .log_interval = 20} };

	auto re_wf = Schrodinger::Spin0{ re_capped_grid, re_potential, re_coupling };
	auto p2 = SplitPropagator<MODE::RE, SplitType, decltype(re_wf)>{ {.dt = 0.3}, std::move(re_wf) };
	p2.run(re_outputs, [=](const WHEN when, const ind step, const uind pass, auto& wf)
		   {
			   if (when == WHEN::AT_START)
			   {
				   if (MPI::region == 0)
					//    wf.load(im_output_name);
					   wf.addUsingCoordinateFunction(
						   [](auto... x) -> cxd
						   {
							//    std::size_t i = 2;
							//    double res = 0.0;
							//    ((res = x * 2.0, true) || ...);
							//    ((i++ == 0 ? (res = x * 2.0, true) : false) || ...);
							//    return sin(((x * 1.9) + ...));
							//    return gaussian(0, 2.0, x...) * cxd { cos(res), sin(res) };
							//    return gaussian(8.0, 4.5, x...) * cxd { cos(((x * 2.0) + ...)), sin(((x * 2.0) + ...)) };
							   return gaussian(8.0, 4.0, x...) * cxd { cos(((x * 2.0) + ...)), sin(((x * 2.0) + ...)) };
						   });

				   logUser("wf loaded manually!");
				   wf.save("at0_");
			   }
			   if (when == WHEN::DURING)
			   {
				   if (step == 1)wf.save("at1_");
				//    if (step == 2)wf.save("at2_");
				   if (step == 50)wf.save("at50_");
				   if (step == 100)wf.save("at100_");
				   if (step == 150)wf.save("at150_");
				   if (step == 200)wf.save("at200_");
				   if (step == 250)wf.save("at250_");
				   if (step == 300)wf.save("at300_");
				//    if (step == 350)wf.save("at350_");
			   }
			   if (when == WHEN::AT_END)
				   wf.save("final");
		   });
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
// 	SinPulse{{.field = 0.1, .omega = 0.06, .ncycles = 4.0, .FWHM_percent = 0.5, .phase_in_pi_units = 0, .delay_in_cycles = 0}}
// };

// using F2 = Field <AXIS::Y, GaussianEnvelope<SinPulse>>;
// F2 f2
// {
// 	GaussianEnvelope<SinPulse>{ {.field = 0.1, .omega = 0.06, .ncycles = 4.0, .FWHM_percent = 0.5, .phase_in_pi_units = 0, .delay_in_cycles = 0}}
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

