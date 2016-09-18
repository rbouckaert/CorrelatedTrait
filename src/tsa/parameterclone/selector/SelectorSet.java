package tsa.parameterclone.selector;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;

@Description("A calculation node that propagates the used values from a vector of parameters. As opposed to Selector, this class aggregates values into some kind of multiset, useful for analyses where only different values matter.")
@Citation("Huelsenbeck, J.P., Larget, B., Alfaro, M.E., 2004. "
		+ "Bayesian Phylogenetic Model Selection Using Reversible Jump Markov Chain Monte Carlo. "
		+ "Mol Biol Evol 21, 1123-1133. doi:10.1093/molbev/msh123")
public class SelectorSet extends Selector {
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			Validate.REQUIRED);
	
	@Override
	public int getDimension() {
		int dim = 0;
		Integer[] sizes = sizesInput.get().getValues();
		for (int i = 0; i<maxIndex; ++i){
			if (sizes[i] > 0) {
				++dim;
			}
		}
		return dim;
	}
	
	@Override
	public double getArrayValue() {
		return -1;
	}

	@Override
	public double getArrayValue(int iDim) {
		int dim = 0;
		Integer[] sizes = sizesInput.get().getValues();
		for (int i = 0; i<maxIndex; ++i){
			if (sizes[i] > 0) {
				++dim;
			}
			if (dim>iDim) {
				int index = groupingsInput.get().getNativeValue(i);
				return parametersInput.get().getValue(index);
			}
		}
		return 0;
	}


}
