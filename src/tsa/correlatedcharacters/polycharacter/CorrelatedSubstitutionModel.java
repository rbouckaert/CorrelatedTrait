/* 
 * Copyright (C) 2015 Gereon Kaiping <gereon.kaiping@soton.ac.uk>
 *
 * This file is part of the BEAST2 package correlatedcharacters.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package tsa.correlatedcharacters.polycharacter;


import java.io.PrintStream;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.evolution.substitutionmodel.Frequencies;
import tsa.parameterclone.selector.Selector;

@Description("Specifies transition probability matrix for a collection of multiple characters."
		+ " At every infinitesimal time step, only one component can change values, so some transition rates are 0, the others arbitrary"
		+ " with restrictions on the rates such that"
		+ " one of the is equal to one and the others are specified relative to"
		+ " this unit rate. Works for any number of states (up to memory use)."
		+ "\nNOTE: While the commonly used DNA substitution models all have the unusual feature that the"
		+ " stationary frequency of the process is built directly into the rate matrix, Pagel&Meade do not"
		+ " assume this for the character substitution model. To use this code correspondingly,"
		+ " it is therefore necessary to set all `frequencies` to equal values.")
@Citation("Pagel, M., Meade, A., 2006." + " Bayesian Analysis of Correlated Evolution of Discrete Characters"
		+ " by Reversible-Jump Markov Chain Monte Carlo."
		+ " The American Naturalist 167, 808--825. doi:10.1086/503444")
public class CorrelatedSubstitutionModel extends ComplexSubstitutionModel implements Loggable {
	// One of these three Inputs is required
	public Input<IntegerParameter> shapeInput = new Input<IntegerParameter>("shape", "component parameter dimensions");
	public Input<CompoundDataType> datatypeInput = new Input<CompoundDataType>("datatype",
			"corresponding compound data type");
	public Input<CompoundAlignment> alignmentInput = new Input<CompoundAlignment>("alignment",
			"corresponding alignment to derive parameter dimensions from");

	protected Integer[] shape;
	
	Function rates;
	protected int nonzeroTransitions = 0;

	public CorrelatedSubstitutionModel() {}
	public CorrelatedSubstitutionModel(IntegerParameter shape, CompoundAlignment characters, Selector rates, Frequencies freqs) {
		initByName("shape", shape,
				"alignment", characters,
				"rates", rates,
				"frequencies", freqs);
	}

	@Override
	public void initAndValidate() {
		frequencies = frequenciesInput.get();

		// One of these three Inputs is required, but there is no way to specify
		// that using Validate.XXX – Therefore we do it manually.
		if (shapeInput.get() == null) {
			CompoundDataType datatype;
			if (alignmentInput.get() == null) {
				if (datatypeInput.get() == null) {
					// FIXME: Construct a ValidateException and use that for
					// cases like this
					throw new IllegalArgumentException("One of shape, datatype, alignment must be specified.");
				} else {
					datatype = datatypeInput.get();
				}
			} else {
				datatype = (CompoundDataType) alignmentInput.get().getDataType();
			}
			shape = datatype.getStateCounts();
		} else {
			shape = shapeInput.get().getValues();
		}

		updateMatrix = true;
		nrOfStates = 1;
		for (int size : shape) {
			nrOfStates *= size;
			nonzeroTransitions += size - 1;
		}

		if (nrOfStates != frequencies.getFreqs().length) {
			throw new RuntimeException("Dimension of input 'frequencies' is " + frequencies.getFreqs().length
					+ " but the " + "shape input gives a total dimension of " + nrOfStates);
		}

		if (ratesInput.get().getDimension() != nrOfStates * nonzeroTransitions) {
			throw new RuntimeException("Dimension of input 'rates' is " + ratesInput.get().getDimension() + " but a "
					+ "rate matrix of dimension " + nrOfStates + "x" + nonzeroTransitions + "="
					+ nrOfStates * nonzeroTransitions + " was " + "expected");
		}

		eigenSystem = createEigenSystem();
//		try {
//			eigenSystem = createEigenSystem();
//			// eigenSystem = new DefaultEigenSystem(m_nStates);
//		} catch (ClassNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (InstantiationException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IllegalAccessException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (InvocationTargetException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

		rateMatrix = new double[nrOfStates][nrOfStates];
		relativeRates = new double[ratesInput.get().getDimension()];
		storedRelativeRates = new double[ratesInput.get().getDimension()];
		rates = ratesInput.get();
	} // initAndValidate

	public Integer[] getShape() {
		return shape.clone();
	}

	/**
	 * sets up rate matrix *
	 */
	protected void setupRateMatrix() {
		// Reset the rate matrix to zero. This is important, because
		// DefaultEigenSystem overwrites it, and sets some zero entries to
		// non-zero.
		rateMatrix = new double[rateMatrix.length][rateMatrix[0].length];

		double[] fFreqs = frequencies.getFreqs();

		rateMatrix[0][1] = this.rates.getArrayValue(0);
		rateMatrix[0][2] = this.rates.getArrayValue(1);
		rateMatrix[1][0] = this.rates.getArrayValue(2);
		rateMatrix[1][3] = this.rates.getArrayValue(3);
		rateMatrix[2][0] = this.rates.getArrayValue(4);
		rateMatrix[2][3] = this.rates.getArrayValue(5);
		rateMatrix[3][1] = this.rates.getArrayValue(6);
		rateMatrix[3][2] = this.rates.getArrayValue(7);

		// bring in frequencies
		for (int i = 0; i < nrOfStates; i++) {
			for (int j = i + 1; j < nrOfStates; j++) {
				rateMatrix[i][j] *= fFreqs[j];
				rateMatrix[j][i] *= fFreqs[i];
			}
		}

		// set the diagonal
		for (int i = 0; i < nrOfStates; i++) {
			double rowsum = 0;
			for (int j = 0; j < nrOfStates; j++) {
				if (j != i) {
					rowsum += rateMatrix[i][j];
				}
			}
			rateMatrix[i][i] = -rowsum;
		}

		// normalise rate matrix to one expected substitution per unit time
		double fSubst = 0.0;
		for (int i = 0; i < nrOfStates; i++) {
			fSubst += -rateMatrix[i][i] * fFreqs[i];
		}
		for (int i = 0; i < nrOfStates; i++) {
			for (int j = 0; j < nrOfStates; j++) {
				rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
			}
		}
		// System.out.println(">" + Arrays.deepToString(rateMatrix));
	} // setupRateMatrix

	public boolean depends(int component, int dependsOn) {
		// Check whether the evolution rates of `component` depend on the state
		// of `dependsOn`.

		// TODO: Currently, this gives the more generic result (“is dependent
		// on”) when `frequencies` are not all equal, but compatible with
		// independent evolution, and the base rates make up for the difference
		// in actual rates introduced by fixing frequencies.
		double[] fFreqs = frequencies.getFreqs();
		double freq0 = fFreqs[0];
		for (double freq : fFreqs) {
			if (freq != freq0) {
				return true;
			}
		}

		Function rates = ratesInput.get();

		int componentMax = 0;
		int componentMin = 0;
		int dependsOnStep = nrOfStates;
		for (int c = 0; c <= component || c <= dependsOn; ++c) {
			if (c <= component) {
				componentMin = componentMax;
				componentMax += shape[c] - 1;
			}
			if (c <= dependsOn) {
				dependsOnStep /= shape[c];
			}
		}

		boolean[] checked = new boolean[nrOfStates];
		for (int from = 0; from < nrOfStates; ++from) {
			if (!checked[from]) {
				checked[from] = true;
				for (int componentTo = componentMin; componentTo < componentMax; ++componentTo) {
					double thisRate = rates.getArrayValue(from * nonzeroTransitions + componentTo);
					for (int other = 1; other < shape[dependsOn]; ++other) {
						int otherIndex = from + dependsOnStep * other;
						checked[otherIndex] = true;
						// System.out.printf("> %d ?= %d\n", from, otherIndex);
						double otherRate = rates.getArrayValue(otherIndex * nonzeroTransitions + componentTo);
						if (thisRate != otherRate) {
							// System.out.printf(" No: %f vs. %f\n", thisRate,
							// otherRate);
							return true;
						}
					}
				}
			}
		}

		return false;
	} // depends


	@Override
	public void init(PrintStream out) {
		out.print("rate_00->01\t");
		out.print("rate_00->10\t");
		out.print("rate_01->00\t");
		out.print("rate_01->11\t");
		out.print("rate_10->00\t");
		out.print("rate_10->11\t");
		out.print("rate_11->01\t");
		out.print("rate_11->10\t");
	}

	@Override
	public void log(int sample, PrintStream out) {
		for (int i = 0; i < 8; i++) {
			out.print(this.rates.getArrayValue(i) + "\t");
		}
	}
	
	@Override
	public void close(PrintStream out) {
	}

} // class GeneralSubstitutionModel
