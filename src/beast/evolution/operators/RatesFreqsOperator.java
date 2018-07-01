package beast.evolution.operators;

import java.text.DecimalFormat;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

@Description("Operator on rates of a asymetric substition model that proposes changes such that frequencies remain the same")
public class RatesFreqsOperator extends Operator {
	
    final public Input<RealParameter> ratesInput =
            new Input<>("rates", "Rate parameter which defines the transition rate matrix. " +
                    "Only the off-diagonal entries need to be specified (diagonal makes row sum to zero in a " +
                    "rate matrix). Entry i specifies the rate from floor(i/(n-1)) to i%(n-1)+delta where " +
                    "n is the number of states and delta=1 if floor(i/(n-1)) <= i%(n-1) and 0 otherwise.", Validate.REQUIRED);

    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    RealParameter rates;
	double scaleFactor = 0.1;
	int [] pair;
	double upper, lower;
	
    @Override
    public void initAndValidate() {
    	rates = ratesInput.get();
    	upper = rates.getUpper();
    	lower = Math.max(rates.getLower(), 0);
    	scaleFactor = scaleFactorInput.get();

		int n = (int) Math.sqrt(rates.getDimension()) + 1;
		int [][] map = new int[n][n];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					map[i][j] = k++;
				}
			}			
		}
		pair = new int[n*(n-1)];
		k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					pair[k++] = map[j][i];
				}
			}			
		}
	}
    
    
    @Override
    public double proposal() {
		int i = Randomizer.nextInt(rates.getDimension());
		int j = pair[i];
		double epsilon = Randomizer.nextDouble() * scaleFactor;
		if (Randomizer.nextBoolean()) {
			// up
			if (rates.getValue(i) + epsilon > upper || rates.getValue(j) + epsilon > upper) {
				return Double.NEGATIVE_INFINITY;
			}
			rates.setValue(i, rates.getValue(i) + epsilon);
			rates.setValue(j, rates.getValue(j) + epsilon);
		} else {
			// down
			double delta = epsilon/(1+epsilon);
			if (rates.getValue(i) - delta < lower || rates.getValue(j) - delta < lower) {
				return Double.NEGATIVE_INFINITY;
			}
			rates.setValue(i, rates.getValue(i) - delta);
			rates.setValue(j, rates.getValue(j) - delta);
		}
		normalise(rates);
		return 0;
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
		
//		sum = 0;
//		for (int i = 0; i < rates.getDimension(); i++) {
//			sum += rates.getValue(i);
//		}
//		assertEquals(sum, rates.getDimension(), 1e-10);
	}

	/**
	 * automatic parameter tuning *
	 */
	@Override
	public void optimize(final double logAlpha) {
	    if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        delta += Math.log(1.0 / scaleFactor - 1.0);
	        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
	    }
	}
	
	public double scaleFactor() {
	    return scaleFactor;
	}
	
	@Override
	public void setCoercableParameterValue(final double value) {
		scaleFactor = Math.max(Math.min(value, upper), lower);
	}
	
	@Override
	public String getPerformanceSuggestion() {
	    final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
	    final double targetProb = getTargetAcceptanceProbability();
	
	    double ratio = prob / targetProb;
	    if (ratio > 2.0) ratio = 2.0;
	    if (ratio < 0.5) ratio = 0.5;
	
	    // new scale factor
	    final double sf = Math.pow(scaleFactor, ratio);
	
	    final DecimalFormat formatter = new DecimalFormat("#.###");
	    if (prob < 0.10) {
	        return "Try setting scaleFactor to about " + formatter.format(sf);
	    } else if (prob > 0.40) {
	        return "Try setting scaleFactor to about " + formatter.format(sf);
	    } else return "";
	}
}
