/* Copyright (C) 2018-2019 Bho Matthiesen, Christoph Hellings
 * 
 * This program is used in the article:
 *
 * Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang
 * Utschick, "Mixed Monotonic Programming for Fast Global Optimization,"
 * submitted to IEEE  Transactions on Signal Processing.
 * 
 * 
 * License:
 * This program is licensed under the GPLv2 license. If you in any way use this
 * code for research that results in publications, please cite our original
 * article listed above.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. */


#include <array>
#include <numeric>
#include <limits>
#include <memory>

#include "ASMAPEL_PA.h"

#include "Gurobi.h"
#include "gurobi_c++.h"

extern "C" {
	#include "mkl.h"
}


constexpr size_t uDIM = 4;
constexpr size_t DIM = (uDIM+1)*(uDIM+1);


using vtype = smapel::PA<DIM>::vtype;

template <size_t Dim>
struct userdata
{
	double N;
	std::array<double, Dim> Pmax;
	std::array<double, Dim> sigma;
	std::array<double, Dim> alpha;
	std::array<std::array<double, Dim>, Dim> beta;

	// computed by preprocess()
	std::array<std::array<double, Dim>, Dim> G;
	std::array<double, Dim> eta;

	constexpr size_t dim() { return Dim; }
	void preprocess(void); // call after setting alpha, beta and sigma
};

template <size_t Dim>
class ASMAPEL : public smapel::PA<Dim>
{
	using typename smapel::PA<Dim>::vtype;
	
	public:
		ASMAPEL(void *ud)
			:
				smapel::PA<Dim>(ud),
				Grb(GurobiEnv::getInstance()),
				Grb_p(Grb.addVars(nullptr, static_cast<userdata<uDIM>*>(ud)->Pmax.data(), nullptr, nullptr, nullptr, uDIM)),
				Grb_t(Grb.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "t")),
				Grb_c(Grb.addConstrs(uDIM)),
				Grb2(GurobiEnv::getInstance()),
				Grb2_beta(Grb2.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, uDIM+1)),
				Grb2_t(Grb2.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "t")),
				Grb2_c(Grb2.addConstrs(uDIM+1))
		{
			/* create gurobi models */
			Grb.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
			Grb2.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
			Grb.set(GRB_IntParam_LogToConsole, 0);
			Grb2.set(GRB_IntParam_LogToConsole, 0);

			// add constraints to Grb
			for (size_t i = 0; i < uDIM; ++i)
			{
				Grb_c[i].set(GRB_CharAttr_Sense, '>');
				Grb.chgCoeff(Grb_c[i], Grb_t, -1.0);
			}

			// add constraints to Grb2
			for (size_t i = 0; i < uDIM+1; ++i)
			{
				Grb2_c[i].set(GRB_CharAttr_Sense, '<');
				Grb2_c[i].set(GRB_DoubleAttr_RHS, 1.0);
				Grb2.chgCoeff(Grb2_c[i], Grb2_beta[i], -1.0);
			}

			GRBLinExpr lhs = 0;
			for (size_t i = 0; i < uDIM+1; ++i)
				lhs += Grb2_beta[i];
			Grb2.addConstr(lhs, '<', 1.0);
		}

		double projTol = 1e-6;

	private:
		typename smapel::PA<Dim>::projection_return projection(const vtype& x, void*) override;

		GRBModel Grb;
		std::unique_ptr<GRBVar[]> Grb_p;
		GRBVar Grb_t;
		std::unique_ptr<GRBConstr[]> Grb_c;

		GRBModel Grb2;
		std::unique_ptr<GRBVar[]> Grb2_beta;
		GRBVar Grb2_t;
		std::unique_ptr<GRBConstr[]> Grb2_c;
};

template <size_t Dim>
void
userdata<Dim>::preprocess(void)
{
	for (size_t i = 0; i < Dim; ++i)
	{
		for (size_t j = 0; j < Dim; ++j)
		{
			if (i == j)
			{
				G[i][j] = 0.0;

				if (beta[i][i] != 0.0)
					throw std::domain_error("beta_ii must be zero");
			}
			else
				G[i][j] = beta[i][j] / alpha[i];
		}

		eta[i] = sigma[i] / alpha[i];
	}
}

template <size_t Dim>
typename smapel::PA<Dim>::projection_return
ASMAPEL<Dim>::projection(const typename ASMAPEL<Dim>::vtype& x, void *user_data)
{
	typename smapel::PA<Dim>::projection_return ret;
	userdata<uDIM>* ud = (userdata<uDIM> *) user_data;
	std::array<double, uDIM+2> lambda;



	for (size_t k = 0; k < uDIM+1; ++k)
	{
		// initialize
		auto p = ud->Pmax;
		double lbd = std::numeric_limits<double>::infinity();

		do
		{
			/* update lambda */
			lbd = std::numeric_limits<double>::infinity();
			for (size_t i = 0; i < ud->dim(); ++i)
			{
				double gamma = ud->alpha[i] * p[i] / std::inner_product(ud->beta[i].begin(), ud->beta[i].end(), p.begin(), ud->sigma[i]);
				double tmp = (1.0 + gamma) / (1.0 + x[i+k*uDIM]);

				lbd = std::min(lbd, tmp);
			}

			/* update p */
			for (size_t i = 0; i < ud->dim(); ++i)
			{
				double tmp = 1 - lbd * (x[i+k*uDIM] + 1);

				Grb_c[i].set(GRB_DoubleAttr_RHS, -tmp*ud->sigma[i]);

				Grb.chgCoeff(Grb_c[i], Grb_p[i], ud->alpha[i]);

				for (size_t j = 0; j < ud->dim(); ++j)
					if (i != j)
						Grb.chgCoeff(Grb_c[i], Grb_p[j], tmp*ud->beta[i][j]);
			}
			Grb.update();
			Grb.optimize();

			for (size_t i = 0; i < ud->dim(); ++i)
				p[i] = Grb_p[i].get(GRB_DoubleAttr_X);

		} while (std::abs(Grb.get(GRB_DoubleAttr_ObjVal)) > projTol);

		lambda[k] = lbd;
		for (size_t i = 0; i < ud->dim(); ++i)
		{
			ret.cand[i+k*uDIM] = ud->alpha[i] * p[i] / std::inner_product(ud->beta[i].begin(), ud->beta[i].end(), p.begin(), ud->sigma[i]);
		}
	}


	/* beta */
	for (size_t k = 0; k < uDIM+1; ++k)
		Grb2.chgCoeff(Grb2_c[k], Grb2_t, 1.0+x[uDIM*(uDIM+1)+k]);

	Grb2.update();
	Grb2.optimize();

	lambda[uDIM+1] = Grb2.get(GRB_DoubleAttr_ObjVal);
	for (size_t k = 0; k < uDIM+1; ++k)
	{
		ret.cand[k+uDIM*(uDIM+1)] = Grb2_beta[k].get(GRB_DoubleAttr_X);
	}

	// compute final lambda
	double lbd = *std::min_element(lambda.begin(), lambda.end());

	// compute projection
	for (size_t i = 0; i < x.size(); ++i)
		ret.ret[i] = lbd*(x[i]+1.0) - 1.0;

	return ret;
}



bool
inH(const vtype& x, void*)
{
	/* check that all are >= 0.0 */
	bool ret = true;
	for (auto it = x.begin(); ret && it != x.end(); ++it)
		ret = ret && *it >= 0.0;
	
	return ret;
}

double
obj(const vtype& x, void*)
{
	double ret = 0;
	
	std::array<double, uDIM> r = {0};
	for (size_t k = 0; k < uDIM+1; ++k)
	{
		for (size_t i = 0; i < uDIM; ++i)
		{
			r[i] += std::max(x[uDIM*(uDIM+1)+k],0.0) * std::log2(1.0 + std::max(x[i+k*uDIM],0.0));
		}
	}
	for (size_t i = 0; i < uDIM; ++i)
	{
		ret += std::log(r[i]);
	}

	return ret;
}

template <size_t Dim>
bool
inG(const vtype& x, void* user_data)
{
	userdata<uDIM>* ud = (userdata<uDIM> *) user_data;

	// TODO how can the following loop be replaced by accumulate?
	//double timesum = std::accumulate(x[uDIM*(uDIM+1)], x[(uDIM+1)*(uDIM+1)-1], 0.0);
	double timesum = 0.0;
	for (size_t k = 0; k < uDIM+1; ++k)
	{
		timesum += std::max(x[k+uDIM*(uDIM+1)],0.0);
	}
	if (timesum > 1.0)
	{
		return false;
	}
	for (size_t k = 0; k < uDIM+1; ++k)
	{
		double G[uDIM][uDIM] __attribute__((aligned(64)));
		double eta[uDIM] __attribute__((aligned(64)));

		for (size_t i = 0; i < uDIM; ++i)
		{
			double gamma = std::max(x[i+k*uDIM],0.0);

			eta[i] = gamma * ud->eta[i];

			for (size_t j = 0; j < uDIM; ++j)
			{
				if (i == j)
					G[i][i] = 1.0;
				else
					G[i][j] = -gamma * ud->G[i][j];
			}
		}

		{
			lapack_int ipiv[uDIM] __attribute__((aligned(64)));

			// WARNING: eta is overwritten, G is overwritten if info != 0
			auto info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, uDIM, 1, &G[0][0], uDIM, ipiv, eta, 1);

			if (info != 0)
				throw std::runtime_error("linear system solver failed with status " + std::to_string(info));
		}

		for (auto&& e : eta)
		{
			e = (abs(e)<1e-14) ? 0 : e;
		}

		for (size_t i = 0; i < uDIM; ++i)
		{
			bool tmp = eta[i] >= 0 && eta[i] <= ud->Pmax[i];

			if (!tmp)
			{
				return false;
			}
		}
	}
	return true;

	/*
	gamma = 2**(x) - 1;
	G_{kj} = gamma_k * beta_{kj} / alpha{k};
	G_{jj} = 0;

	eta_k = gamma_k * sigma[k] / alpha[k];

	Solve (I-G) *p = eta for p
	
	return all(0 <= p <= pmax)
	*/
}

int
main(void)
{
/*
	const double Pmax = 1;
	const double N = 1e-2;

	const double alpha[] = { 5.0091360e-0, 2.4832193e+00, 2.6642494e+00, 1.9597458e+00, 2.6108010e-01, 1.3730196e+00, 3.9751133e-01, 1.1865231e-01, 6.8428694e+00, 4.3530133e+00, 1.2695716e+00, 5.4620505e+00, 8.1000105e-01, 2.5237667e-01, 4.3768792e+00, 1.1888159e+00, 9.0684972e-02, 1.4780348e+00, 2.2194814e+00, 2.4494290e+00 };
	const double beta[][20] = {
		{ 3.3516109e-01, 7.9828486e-01, 1.0228757e+00, 3.5630695e-01, 5.0516285e+00, 4.0510208e-01, 1.2883211e-01, 1.0126427e+00, 1.3000776e-01, 9.9310068e-01, 1.0595030e+00, 2.5490221e+00, 7.0220074e-01, 1.8244044e+00, 6.0104347e-01, 1.1408571e+00, 1.9786585e-01, 3.4937312e+00, 1.8277113e-01, 2.0552713e-01 },
		{ 1.7795534e-01, 1.6365220e+00, 8.9789916e-02, 7.9793004e-01, 3.0525077e+00, 7.9354112e-02, 4.5552998e-01, 2.7003983e-01, 1.5274112e+00, 1.3211601e+00, 1.1239266e+00, 4.1683973e+00, 3.3579536e-01, 4.5727826e-01, 1.7165907e-01, 6.3784642e-02, 8.8094221e-02, 7.4063286e-01, 1.3093111e+00, 2.6932336e-01 },
		{ 6.0391173e-01, 7.3569861e-01, 2.8252909e-02, 2.1882761e-02, 1.3575220e+00, 2.0614186e+00, 2.6341307e+00, 3.7052207e-01, 4.5076661e-01, 9.6808104e-01, 1.5499129e+00, 3.0410614e-01, 2.1540984e-01, 6.3280160e-01, 7.0377232e-02, 1.2448899e+00, 1.0981851e+00, 4.8820963e-01, 1.1204575e+00, 2.9384772e-03 },
		{ 1.2793026e-01, 6.2550964e-01, 1.4376196e+00, 1.9053533e-01, 2.2446574e+00, 1.0033350e-01, 3.6425315e+00, 1.3986873e+00, 4.8289251e-02, 1.1688305e-01, 4.1958383e-02, 9.1200694e-01, 9.3656497e-01, 2.5489586e-01, 3.0849290e-01, 3.9350019e-01, 4.9384952e-01, 8.9714932e-01, 6.4545232e-01, 1.7061307e+00 },
		{ 6.5668148e-01, 2.3490561e+00, 8.7151218e-01, 5.2440691e-02, 1.1202852e-01, 8.5846895e-01, 1.4636498e+00, 4.4165312e-01, 1.4564059e+00, 1.9985075e+00, 1.2873553e+00, 2.9167385e-01, 8.3286798e-01, 3.8105600e-01, 9.8959494e-02, 9.3236808e-01, 8.9347466e-01, 6.6062458e-01, 1.7325593e-01, 2.6620000e-01 },
		{ 1.1499842e+00, 2.6551308e-01, 2.6223660e-01, 2.7697495e-01, 3.7749186e-02, 6.7745350e-01, 2.5785098e-01, 8.7506576e-01, 2.6264241e-01, 1.8554610e-01, 1.2307743e+00, 1.4452999e+00, 1.1251780e-01, 4.7693628e-01, 2.1252393e-01, 1.2428201e+00, 5.1850858e-01, 5.6526761e+00, 9.0899080e-02, 1.5369746e-01 },
		{ 1.4122834e-02, 1.2883938e+00, 4.0237261e-01, 1.2520063e-01, 3.0834071e+00, 6.9132011e-01, 8.5500076e-01, 7.5961717e-01, 3.4170709e-02, 7.6208894e-01, 3.1331044e-01, 1.2074809e-01, 2.5124150e-01, 2.5540292e-01, 8.3422210e-01, 4.7952879e-01, 1.8898094e-01, 7.0714292e-01, 7.1923509e-01, 2.0611055e+00 },
		{ 1.5352557e+00, 3.5604144e-01, 5.8092338e-02, 4.1579478e+00, 4.1001712e-02, 5.1705751e-01, 6.2365429e-01, 7.9473044e-01, 8.6627545e-01, 3.9056632e-02, 1.2863244e+00, 9.8807134e-01, 2.9603841e-01, 6.3254917e-01, 7.3956164e-01, 1.3534253e-01, 8.8662307e-01, 6.0877074e-01, 1.0820750e-01, 2.6168716e-01 },
		{ 6.8795227e-01, 1.0072855e+00, 4.9614160e-02, 2.1458335e+00, 3.7476855e-01, 2.5493124e+00, 1.0178964e+00, 8.8551854e-01, 2.4398088e-02, 2.0371166e-01, 9.1707630e-03, 3.0405780e-01, 1.1643111e+00, 1.5157698e+00, 2.8114443e+00, 4.6637007e-01, 5.5375367e-02, 7.3355588e-01, 9.7531860e-01, 1.8976760e+00 },
		{ 7.8895947e-01, 2.8074705e+00, 7.7639511e-01, 1.7380545e+00, 2.1327782e+00, 5.0133777e-01, 1.6004050e+00, 2.3482287e-03, 9.9758497e-01, 6.0483270e-01, 7.2357268e-02, 4.0647120e-01, 1.3493438e+00, 7.4638903e-01, 1.0137375e+00, 1.3329636e+00, 1.3154425e-01, 6.8774837e-01, 2.8682197e-01, 4.0759441e-01 },
		{ 3.7491947e-01, 2.5907370e-01, 8.4915958e-01, 2.9368584e+00, 9.2814923e-02, 2.3121518e-01, 1.4163872e+00, 2.5697495e-01, 6.7622173e-01, 7.8871741e-02, 1.3686664e+00, 2.6303512e-01, 3.3427597e-03, 1.4053089e+00, 2.1282719e-01, 1.4538582e+00, 1.7915503e+00, 1.1704038e+00, 1.0535411e+00, 6.1475861e-01 },
		{ 2.0471512e+00, 6.1232019e-01, 4.3266983e-02, 9.5015459e-01, 3.4160055e-01, 8.3707905e-02, 3.7769755e+00, 2.4592014e-01, 1.4572429e+00, 9.6227408e-01, 1.4347772e-03, 6.7662882e-01, 2.7760011e+00, 1.0478128e-02, 1.4550473e+00, 2.1535710e+00, 2.2775271e-01, 3.0585588e-01, 1.7248772e+00, 4.6169419e-02 },
		{ 1.2193313e+00, 3.3762265e-01, 4.4148727e-01, 1.0114714e+00, 1.4287395e+00, 4.8428513e-01, 1.0510682e+00, 4.1931423e-01, 6.8585845e-01, 1.3761764e+00, 4.8858511e-01, 1.0558564e+00, 5.4342893e-01, 9.5397796e-02, 7.6386119e-02, 3.1470282e-01, 4.0060759e-01, 4.6502415e-01, 2.3719088e-01, 2.1867175e-01 },
		{ 8.0752700e-01, 4.0203890e-01, 3.3466429e+00, 5.8894953e-01, 8.7638641e-01, 3.6284435e-01, 1.5803210e+00, 6.6975347e-03, 2.7761073e-01, 3.3502944e-01, 1.6016509e+00, 4.3010737e-01, 3.8031975e-01, 5.7551957e-01, 6.2435501e+00, 2.1196391e+00, 5.5187268e-03, 5.7598826e-01, 1.1338238e+00, 3.8086809e-01 },
		{ 1.5297003e-03, 4.5565910e-01, 1.7185472e+00, 4.8559463e-01, 8.9843745e-02, 1.1728606e+00, 8.8651462e-02, 8.0661829e-02, 1.2092959e+00, 1.4547281e-01, 1.6203188e-01, 3.1800086e+00, 6.8603790e-01, 3.1810666e+00, 5.9095770e-01, 2.4445121e+00, 3.7090258e-01, 1.4020798e+00, 2.2514341e-01, 2.3926240e+00 },
		{ 1.8002585e+00, 1.0683794e+00, 2.0271754e-01, 1.5421566e-01, 9.0220934e-01, 1.5018342e-01, 8.6723971e-01, 7.3189125e-01, 7.1927123e-01, 4.9274712e-01, 1.3520450e+00, 1.0113979e+00, 9.6198531e-01, 4.3116641e-01, 1.1745327e+00, 3.5834901e-01, 1.6112471e+00, 1.0247893e-01, 4.4888611e-01, 1.0267920e+00 },
		{ 4.9236785e-01, 1.4986235e+00, 2.7503122e-01, 1.1078127e+00, 1.4220507e+00, 7.2821873e-01, 2.4141956e+00, 1.1820753e+00, 3.4526467e-01, 7.8210863e-04, 2.0614981e+00, 1.5427963e+00, 3.1390116e+00, 9.0374619e-01, 2.8859152e+00, 1.1570203e-01, 5.6266804e-02, 1.8752037e+00, 4.5347153e-02, 1.8352951e+00 },
		{ 1.0009066e-01, 4.0763867e-01, 3.3993276e+00, 2.2047369e+00, 1.5216873e+00, 1.0494945e+00, 7.2450978e-02, 2.6072490e-01, 6.0250628e-02, 6.3859430e+00, 2.6208693e+00, 4.2905488e-01, 2.0000303e-01, 7.7201660e-02, 6.7079317e-01, 5.5141094e-01, 6.3202603e-02, 1.4738896e+00, 5.8642963e-01, 6.9726005e-01 },
		{ 5.1848753e-01, 2.3168775e-01, 5.1083491e-01, 1.6006393e+00, 1.4905275e-01, 4.0850873e-02, 3.6793287e-01, 2.5520994e-01, 1.5574769e+00, 2.3685738e+00, 2.5109673e+00, 8.2996925e-01, 2.5302645e+00, 7.3179681e-01, 3.1701467e+00, 9.6586909e-01, 2.4040233e+00, 2.2575332e-01, 6.0648980e-02, 1.5703631e+00 },
		{ 1.0995063e+00, 2.0157725e-01, 1.9863521e+00, 1.4773864e+00, 5.4565338e-02, 6.4464356e-01, 4.2191446e-01, 8.8059381e-01, 2.1095148e+00, 1.4644133e+00, 7.1697435e-01, 7.0410651e-02, 3.0879528e+00, 7.2283272e-01, 1.9877261e+00, 2.3237078e+00, 8.5798775e-01, 2.5310085e+00, 6.3721946e-01, 1.7257242e+00 }
	};
*/

/*
	const double Pmax = 10;
	const double N = 1;
	const double alpha[] = { 4.1250, 0.9870 };
	const double beta[][20] = {
		{ 0, 2.1803 },
		{ 0.5300, 0 }
	};
*/
	
	// Scenario from ASMAPEL paper
	const double Pmax = 1;
	const double N = 0.1e-3;
	const double alpha[] = { 6.2500e-02, 2.4414e-04, 7.7160e-04, 3.9062e-03 };
	const double beta[][20] = {
		{ 0.0, 2.0475e-05, 9.2456e-05, 1.2625e-04 },
		{ 5.9488e-04, 0.0, 1.2625e-04, 2.6031e-05 },
		{ 9.2456e-05, 6.4000e-05, 0.0, 3.5013e-05 },
		{ 2.6874e-04, 9.5260e-06, 1.1891e-03, 0.0 }
	};

	userdata<uDIM> ud;

	for (size_t i = 0; i < ud.dim(); ++i)
	{
		ud.alpha[i] = alpha[i];
		for (size_t j = 0; j < ud.dim(); ++j)
		{
			if (j == i)
				ud.beta[i][j] = 0.0; // beta_i,i is not supported
			else
				ud.beta[i][j] = beta[i][j];
		}
		ud.sigma[i] = N;
		ud.Pmax[i] = Pmax;
	}

	ASMAPEL<DIM> pa(&ud);

	for (size_t k = 0; k < uDIM+1; ++k)
	{	for (size_t i = 0; i < uDIM; ++i)
		{
			double tmp = alpha[i] * Pmax / N;
			pa.setUB(i+k*uDIM, tmp);
			std::cout << tmp << " : ";
		}
		std::cout << " | ";
		pa.setUB(k+uDIM*(uDIM+1), 1);
	}
	std::cout << std::endl;

	pa.setPrecision(1e-2);
	pa.setObjective(obj);
	pa.setInH(inH);
	pa.setInG(inG<DIM>);
	//pa.setShift(1); // do not shift since custom projection accounts for the shift
	pa.output_every = 1;

	ud.preprocess();
	pa.optimize();
}
