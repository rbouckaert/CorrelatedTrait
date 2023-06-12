package beast.evolution.substitutionmodel;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;

@Description("substitution model that allows no change of state, so a state (distribution) at "
		+ "the top of a branch will be exactly the same at the bottom of the branch.")
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
