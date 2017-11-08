package beast.evolution.substitutionmodel;

import java.util.Arrays;

import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

public class NoChangeSubstitutionModel extends SubstitutionModel.Base {

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		nrOfStates = frequencies.getFreqs().length;
	}
	
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		Arrays.fill(matrix, 0.0);
		double [] freqs = frequencies.getFreqs();
		for (int i = 0; i < nrOfStates; i++) {
			matrix[i * nrOfStates + i] = freqs[i] > 0 ? 1.0 : 0.0;
		}
	}

	
	@Override
	public double[] getRateMatrix(Node node) {
		throw new RuntimeException("BEAGLE not supported: run with java only");
	}
	
	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		throw new RuntimeException("BEAGLE not supported: run with java only");
		//return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return false;
	}

}
