package beast.evolution.substitutionmodel;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;

@Description("Loggable/Function returning stationary frequencies of a substitution model. "
		+ "The model must have the getTransitionProbabilities() method implemented.")
public class StationaryFrequenciesLogger extends BEASTObject implements Function, Loggable {
	public Input<SubstitutionModel> substModelInput = new Input<>("substModel", "substitution model for which to calculate stationary frequencies", Validate.REQUIRED);

	SubstitutionModel substModel;
	
	@Override
	public void initAndValidate() {
		substModel = substModelInput.get();
	}

	@Override
	public void init(PrintStream out) {
		int dim = substModel.getStateCount();
		for (int i = 0; i < dim; i++) {
			out.append("statFreqs" + getID() + "." + (i+1) + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		double [] freqs = getStationaryFreqs();
		for (int i = 0; i < freqs.length; i++) {
			out.append(freqs[i] + "\t");
		}
	}
	
	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

	@Override
	public int getDimension() {
		return substModel.getStateCount();
	}

	@Override
	public double getArrayValue() {
		return getArrayValue(0);
	}

	@Override
	public double getArrayValue(int dim) {
		double [] freqs = getStationaryFreqs();
		return freqs[dim];
	}

	private double[] getStationaryFreqs() {
		int n = substModel.getStateCount();
		double [] matrix = new double[n*n];
		double [] freqs = new double[n];
		substModel.getTransitionProbabilities(null, 1000, 0, 1.0, matrix);
		System.arraycopy(matrix, 0, freqs, 0, n);
		return freqs;
	}

}
