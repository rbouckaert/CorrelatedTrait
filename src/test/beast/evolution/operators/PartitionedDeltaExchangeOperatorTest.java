package test.beast.evolution.operators;

import java.util.Arrays;

import org.junit.Test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SVSGeneralSubstitutionModel;
import junit.framework.TestCase;

public class PartitionedDeltaExchangeOperatorTest extends TestCase {

	
	@Test
	public void testQ() {
		BooleanParameter mask = new BooleanParameter();
		mask.initByName("dimension", 20, "value", true);

		RealParameter rates = new RealParameter();
		rates.initByName("value", 
				       "1.0 1.0 1.0 1.0"
				+ " 2     1.0 1.0 1.0 "
				+ " 0 1.0     1.0 1.0 "
				+ " 1.0 1.0 1.0     1.0 "
				+ " 1.0 1.0 1.0 1.0 "
				+ "");
		
		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", "0.3 0.15 0.15 0.2 0.2");

		SVSGeneralSubstitutionModel Q = new SVSGeneralSubstitutionModel();
		Q.initByName("frequencies", freqs, "rateIndicator", mask, "rates", rates, "symmetric", false);
		
		double [] matrix = new double[36];
		Q.getTransitionProbabilities(null, 100.0, 0.0, 1.0, matrix);
		System.out.println(Arrays.toString(matrix));
		
		double [] d = new double[5];
		for (int i = 0; i < 5; i++) {
			System.arraycopy(matrix, i * 5, d, 0, 5);
			System.out.println(Arrays.toString(d));			
		}
	}
}
