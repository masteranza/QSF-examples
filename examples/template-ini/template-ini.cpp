#include "QSF.h"

constexpr auto opt = OPTIMS::NONE;
constexpr auto order = 1;
using VTV = Split3Base<REP::X, REP::P, REP::X>;
using SplitType = MultiProductSplit<VTV, order>;

using namespace QSF;
using namespace QSF::IO;
int main(int argc, char* argv[])
{

	QSF::init(argc, argv, project_dir / results_dir);
	using grid_t = CartesianGrid<3_D>;
	std::string im_output_name = "groundstate";


	using im_wf_t = Schrodinger::Spin0<grid_t, CoulombInteraction>;
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

	if constexpr MODE_FILTER_OPT(MODE::IM)
	{
		auto p1 = QSF::SplitPropagator<MODE::IM, SplitType, im_wf_t>{};
		p1.run<im_outputs_t>(
			[&](const WHEN when, const ind step, const uind pass, auto& wf)
			{
				if (when == WHEN::AT_START)
				{
					wf.addUsingCoordinateFunction(
						[](auto... x) -> cxd
						{
							return cxd{ gaussian(0.0, 2.0, x...), 0 };
						});
						
					logUser("wf loaded manually!");
				}
				if (when == WHEN::AT_END)
					wf.save(im_output_name);
			});
	}


	if constexpr MODE_FILTER_OPT(MODE::RE)
	{
		using F1 = Field<AXIS::X, GaussianEnvelope<SinPulse>>;
		using re_outputs_t = BufferedBinaryOutputs <
			VALUE<Step, Time>
			, VALUE<F1>
			, AVG<Identity>
			, AVG<PotentialEnergy>
		// , PROJ<EIGENSTATES, Identity>
			, AVG<DERIVATIVE<0, PotentialEnergy>>
			// , FLUX<BOX<3>>
			, VALUE<ETA> /* Estimated end of computation */
		>;
		using coupling_t = DipoleCoupling<VelocityGauge, F1>;
		using re_wf_t = Schrodinger::Spin0<CAP<grid_t>, CoulombInteraction, coupling_t>;

		auto p2 = SplitPropagator<MODE::RE, SplitType, re_wf_t>{};
		p2.run<re_outputs_t>(
			[=](const WHEN when, const ind step, const uind pass, auto& wf)
			{
				if (when == WHEN::AT_START)
					wf.load(im_output_name);

				if (when == WHEN::DURING && step % 50 == 0)
					wf.snapshot("step_" + std::to_string(step));

				if (when == WHEN::AT_END)
					wf.save("final");
			});
	}

	finalize();
}
