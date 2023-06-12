package tsa.parameterclone.helpers;


import org.apache.commons.math.distribution.Distribution;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.distribution.ParametricDistribution;

@Citation("Huelsenbeck, J.P., Larget, B., Alfaro, M.E., 2004. "
		+ "Bayesian Phylogenetic Model Selection Using Reversible Jump Markov Chain Monte Carlo. "
		+ "Mol Biol Evol 21, 1123-1133. doi:10.1093/molbev/msh123")
@Description("Rescaled flat Dirichlet distribution.  p(x_1,...,x_K) = Gamma(K) Prod_{i=1}^k n_i/6^K")
public class RescaledDirichlet extends ParametricDistribution {
	public Input<IntegerParameter> nInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each x");

	public RescaledDirichlet() {}
	public RescaledDirichlet(IntegerParameter sizes) {
		initByName("sizes", sizes);
	}

	@Override
	public double calcLogP(Function pX) {
		double fLogP = 0;
		int K = 0;
		int N = 0;
		for (int i = 0; i < pX.getDimension(); i++) {
			int n_i = nInput.get().getNativeValue(i);
			if (n_i > 0) {
				fLogP += Math.log(n_i);
				K += 1;
			}
			N += n_i;
		}
		fLogP += org.apache.commons.math.special.Gamma.logGamma(K);
		fLogP -= K * Math.log(N);
		return fLogP;
	}

	@Override
	public void initAndValidate() {
	}

	@Override
	public Distribution getDistribution() {
		// TODO Auto-generated method stub
		return null;
	}

}
