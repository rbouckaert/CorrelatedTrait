/* 
 * Copyright (C) 2015 Gereon Kaiping <gereon.kaiping@soton.ac.uk>
 *
 * This file is part of the BEAST2 package parameterclone.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
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

package tsa.parameterclone.selector;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

@Description("A calculation node that propagates parameters from a vector of parameters")
// Needed for the reversible-jump Markov chain described in
@Citation("Huelsenbeck, J.P., Larget, B., Alfaro, M.E., 2004. "
		+ "Bayesian Phylogenetic Model Selection Using Reversible Jump Markov Chain Monte Carlo. "
		+ "Mol Biol Evol 21, 1123-1133. doi:10.1093/molbev/msh123")
public class Selector extends CalculationNode implements Loggable, Function {
	// Input objects
	final public Input<IntegerParameter> entryInput = new Input<IntegerParameter>(
			"entry",
			"The index of the parameter vector that this object propagates",
			(IntegerParameter) null);
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new RealParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", new IntegerParameter(),
			Validate.REQUIRED);
	public Input<Function> minusOneIsSpecialInput = new Input<Function>(
			"minusOne", "Deliver this constant value for index -1"); 

	// Member objects
	protected IntegerParameter entries;
	protected Integer maxIndex;


	public Selector() {		
	}
	
	public Selector(IntegerParameter indices, RealParameter parameters, IntegerParameter groupings) {
		initByName("groupings", indices, "parameters", parameters, "entry", groupings);
	}

	@Override
	public void initAndValidate() {
		try {
			maxIndex = parametersInput.get().getDimension();
		} catch (NullPointerException e) {
			throw new IllegalArgumentException("parameters was left unspecified");
		}
		if (entryInput.get() == null) {
			Integer[] build_entries = new Integer[maxIndex];
			for (int i=0; i<maxIndex; ++i) {
				build_entries[i] = i;
			}
			entries = new IntegerParameter(build_entries);
		} else {
			entries = entryInput.get();
			for (Integer entry : entries.getValues()) {
				if (entry > groupingsInput.get().getDimension()) {
					throw new RuntimeException(
							"entries must be valid index of groupings");
				}
			}
		}
		// Array-like RealParameters do not implement java.lang.iterable, so we
		// must do the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new RuntimeException(
						"All entries in groupings must be valid indices of parameters");
			}
		}
		// value = parametersInput[groupingsInput[entry]]
	}

	/**
	 * Function interface implementation follows *
	 */

	@Override
	public int getDimension() {
		return entries.getDimension();
	}

	@Override
	public double getArrayValue() {
		int index = groupingsInput.get().getNativeValue(
				entries.getValue());
		if (index == -1) {
			return minusOneIsSpecialInput.get().getArrayValue();
		} else {
			return parametersInput.get().getValue(index);
		}
	}

	@Override
	public double getArrayValue(int iDim) {
		int index = groupingsInput.get().getNativeValue(
				entries.getNativeValue(iDim));
		if (index == -1) {
			return minusOneIsSpecialInput.get().getArrayValue();
		} else {
			return parametersInput.get().getValue(index);
		}
	}

	/**
	 * Loggable interface implementation follows *
	 */

	@Override
	public void init(final PrintStream out) {
		for (int i = 0; i < getDimension(); ++i) {
			out.print(getID() + "" + i + "\t");
		}
	}

	@Override
	public void log(final long nSample, final PrintStream out) {
		for (int i = 0; i < getDimension(); ++i) {
			out.print(getArrayValue(i) + "\t");
		}
	}

	@Override
	public void close(final PrintStream out) {
		// nothing to do
	}

}
