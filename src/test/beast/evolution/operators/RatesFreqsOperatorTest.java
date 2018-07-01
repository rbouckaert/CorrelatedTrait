package test.beast.evolution.operators;

import java.text.DecimalFormat;
import java.util.Arrays;

import org.junit.Test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.RatesFreqsOperator;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SVSGeneralSubstitutionModel;
import junit.framework.TestCase;

public class RatesFreqsOperatorTest extends TestCase {

	
	@Test
	public void testQ() {
		BooleanParameter mask = new BooleanParameter();
		mask.initByName("dimension", 20, "value", true);

		RealParameter rates = new RealParameter();
		rates.initByName("value", 
				       "1.0 1.0 1.0 1.0"
				+ " 1.0     1.0 1.0 1.0 "
				+ " 1.0 1.0     1.0 1.0 "
				+ " 1.0 1.0 1.0     1.0 "
				+ " 1.0 1.0 1.0 1.0 "
				+ "");
		
		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", "0.3 0.15 0.15 0.2 0.2");

		SVSGeneralSubstitutionModel Q = new SVSGeneralSubstitutionModel();
		Q.initByName(
//				"eigenSystem","beast.evolution.substitutionmodel.RobustEigenSystem",
				"frequencies", freqs, "rateIndicator", mask, "rates", rates, "symmetric", false);
		
		double [] matrix = new double[36];
		Q.getTransitionProbabilities(null, 100.0, 0.0, 1.0, matrix);
		printQ(matrix);
		
		double [] f = freqs.getFreqs();
		
		int n = 5;
		int [][] map = new int[n][n];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					map[i][j] = k++;
				}
			}			
		}
		int [] pair = new int[n*(n-1)];
		k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					pair[k++] = map[j][i];
				}
			}			
		}
		
		
		
		
		DecimalFormat formatter = new DecimalFormat("#.######");
		for (int i = 0; i < 20; i++) {
			double epsilon = 1e-1;
			double i0 = rates.getValue(i);
			double p0 = rates.getValue(pair[i]);
			rates.setValue(i, i0 + epsilon);
			rates.setValue(pair[i], p0 + epsilon);
			normalise(rates);
			double [] gradient = getGradient(matrix, Q, f);
			
			// should be all zeros
			System.out.print(i+" " + pair[i] +" [");
			for (k = 0; k < 5; k++) {
				System.out.print(formatter.format(gradient[k])+ ", ");
			}
			System.out.println("]");
		}

		
		RatesFreqsOperator operator = new RatesFreqsOperator();
		operator.initByName("rates", rates, "weight", 1.0);
		
		int N = 100000;
		for (int i = 0; i < N; i++) {
			operator.proposal();
		}
		double [] gradient = getGradient(matrix, Q, f);
		// should be all zeros
		System.out.print("diff in frequencies after " + N +" proposals = [");
		for (k = 0; k < 5; k++) {
			System.out.print(formatter.format(gradient[k])+ ", ");
			assertEquals(0.0, gradient[k], 1e-10);
		}
		System.out.println("]");
		
	}
	

	private void normalise(RealParameter rates) {
		double sum = 0;
		for (int i = 0; i < rates.getDimension(); i++) {
			sum += rates.getValue(i);
		}
		sum /= rates.getDimension();
		for (int i = 0; i < rates.getDimension(); i++) {
			rates.setValue(i, rates.getValue(i)/sum);
		}
		
		sum = 0;
		for (int i = 0; i < rates.getDimension(); i++) {
			sum += rates.getValue(i);
		}
		assertEquals(sum, rates.getDimension(), 1e-10);
	}

	private double [] getGradient(double [] matrix, SVSGeneralSubstitutionModel Q, double [] f) {
		Q.makeDirty();
		Q.getTransitionProbabilities(null, 100.0, 0.0, 1.0, matrix);
		double [] gradient = new double[5];
		for (int k = 0; k < 5; k++) {
			gradient[k] = f[k] - matrix[k]; 
		}
		Q.makeDirty();
		return gradient;
	}

	
	private void printQ(double [] matrix) {
		double [] d = new double[5];
		for (int i = 0; i < 5; i++) {
			System.arraycopy(matrix, i * 5, d, 0, 5);
			System.out.println(Arrays.toString(d));			
		}
	}
}
