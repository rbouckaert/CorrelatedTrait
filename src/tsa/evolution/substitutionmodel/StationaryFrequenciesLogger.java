package tsa.evolution.substitutionmodel;


import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.core.Input.Validate;

@Description("Loggable/Function returning stationary frequencies of a substitution model. "
		+ "The model must have the getTransitionProbabilities() method implemented.")
public class StationaryFrequenciesLogger extends CalculationNode implements Function, Loggable {
	public Input<SubstitutionModel> substModelInput = new Input<>("substModel", "substitution model for which to calculate stationary frequencies", Validate.REQUIRED);
	public Input<Integer> dimInput = new Input<>("dim","dimension to be logged. If not specified, all dimensions are logged, otherwise only the first dim dimensions are available", -1);

	SubstitutionModel substModel;
	int dim;
	
	@Override
	public void initAndValidate() {
		substModel = substModelInput.get();
		dim = substModel.getStateCount();
		if (dimInput.get() > 0) {
			dim = dimInput.get();
			if (dim > substModel.getStateCount()) {
				throw new IllegalArgumentException("dim-input cannot be larger than state count (=" + substModel.getStateCount() + ")");
			}
		}
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < dim; i++) {
			out.append("statFreqs" + getID() + "." + (i+1) + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		double [] freqs = getStationaryFreqs();
		for (int i = 0; i < dim; i++) {
			out.append(freqs[i] + "\t");
		}
	}
	
	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

	@Override
	public int getDimension() {
		return dim;
	}

	@Override
	public double getArrayValue() {
		return getArrayValue(0);
	}

	@Override
	public double getArrayValue(int dim) {
		if (dim >= this.dim) {
			throw new IllegalArgumentException("dim-input should be larger entry (=" + dim + ")");			
		}
		double [] freqs = getStationaryFreqs();
		return freqs[dim];
	}

	private double[] getStationaryFreqs() {
		int n = substModel.getStateCount();
		double [] matrix = new double[n*n];
		double [] freqs = new double[dim];
		substModel.getTransitionProbabilities(null, 1000, 0, 1.0, matrix);
		System.arraycopy(matrix, 0, freqs, 0, dim);
		return freqs;
	}

}
