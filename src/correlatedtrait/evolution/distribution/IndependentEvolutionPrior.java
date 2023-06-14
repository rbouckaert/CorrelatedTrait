package correlatedtrait.evolution.distribution;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import correlatedtrait.evolution.substitutionmodel.CorrelatedSubstitutionModel;

@Description("A distribution on correlated substitution models,"
		+ " assuming a prior probability distribution on independent vs." + " dependent evolution models.")
public class IndependentEvolutionPrior extends Distribution {
	public Input<CorrelatedSubstitutionModel> csmInput = new Input<CorrelatedSubstitutionModel>("model",
			"The CorrelatedSubstitutionModel this prior is conditioning");
	public Input<RealParameter> pIndependentInput = new Input<RealParameter>("pIndependent",
			"prior probability for two characters to be independent");

	private CorrelatedSubstitutionModel csm;
	private double logPDependent;
	
	@Override
    public double calculateLogP() {
		if (isDirtyCalculation()) {
			calcLogP();
		}
		return logP;
	}
	
	private void calcLogP() {
		double pIndependent = pIndependentInput.get().getValue();
		logPDependent = Math.log(1-pIndependent) - Math.log(pIndependent);
		
		logP = 0;
		csm = csmInput.get();
		int components = csm.getShape().length;
		for (int c1 = 0; c1 < components; ++c1) {
			for (int c2 = 0; c2<components; ++c2) {
				if (csm.depends(c1, c2)) {
					logP += logPDependent;
				}
				if (csm.depends(c2, c1)) {
					logP += logPDependent;
				}
			}
		}
	}
	
	@Override public List<String> getArguments() { return null; }
	@Override public List<String> getConditions() { return null; }
	@Override public void sample(State state, Random random) { }

}
