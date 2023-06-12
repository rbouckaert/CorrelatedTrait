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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.datatype.DataType;

/**
 * @author Gereon Kaiping <anaphory@yahoo.de>
 * @since 2015-04-24
 */
@Description("Compound datatype. Represents tuples of values from other data types.")
public class CompoundDataType extends DataType.Base {
	public Input<IntegerParameter> reserveAmbiguities = new Input<IntegerParameter>("ambiguities",
			"Number of additional entries to reserve internally for ambiguities",
			new IntegerParameter(new Integer[] {1}),
			Validate.OPTIONAL);
	public Input<List<DataType>> componentsInput = new Input<List<DataType>>("components",
			"Component data types for this compound", new ArrayList<DataType>(), Validate.REQUIRED);
	public Input<IntegerParameter> ambiguitiesSizesInput = new Input<IntegerParameter>(
			"componentSizesIncludingAmbiguities",
			"Number of different specifications of each component – Inferred otherwise."
					+ " This is useful when sub-states include ambiguities. It is never reported externally.",
			(IntegerParameter) null);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"componentSizes",
			"Number of different specifications of each component – Inferred otherwise.",
			(IntegerParameter) null);
	public Input<String> splitInput = new Input<String>("split", "How to split alignment values into separate traits",
			";;");

	protected List<DataType> components;
	protected Integer[] stateCountsIncludingAmbiguities;
	protected Integer[] stateCountsExcludingAmbiguities;
	protected int stateCount = 1;

	public CompoundDataType(List<DataType> inputs, Integer[] sizes, Integer[] sizesWithAmbiguities) {
		super();
		initAndValidate(inputs, new IntegerParameter(sizes), new IntegerParameter(sizesWithAmbiguities));
	}

	public CompoundDataType() {
		super();
	}

	@Override
	public void initAndValidate() {
		initAndValidate(componentsInput.get(), ambiguitiesSizesInput.get(), sizesInput.get());
	}

	private void initAndValidate(List<DataType> components_, IntegerParameter ambiguitiesSizes, IntegerParameter sizes) {
		components = components_;
		if (sizes == null) {
			// We have to rely on the components to tell us their sizes
			stateCountsExcludingAmbiguities = new Integer[components.size()];
			int n = 0;
			for (DataType t : components) {
				int size = t.getStateCount();
				if (size > 0) {
					stateCount *= size;
					stateCountsExcludingAmbiguities[n] = size;
					++n;
				} else {
					throw new IllegalArgumentException("Can only compound finite DataTypes");
				}
			}
		} else {
			// Rely on sizes.
			stateCountsExcludingAmbiguities = sizes.getValues();
			if (components_.size() == 1) {
				// All components are of the same type, yay!
				// Yes, do start this iteration at i=1, because i=0 is already
				// filled.
				DataType component = components.get(0);
				components = new ArrayList<DataType>();
				for (int i = 0; i < sizes.getDimension(); ++i) {
					stateCount *= stateCountsExcludingAmbiguities[i];
				}
			} else {
				if (sizes.getDimension() == components_.size()) {
					// We know the data types, and we don't trust them, getting
					// our data from sizes instead.
					for (int i = 0; i < sizes.getDimension(); ++i) {
						stateCount *= stateCountsExcludingAmbiguities[i];
					}
				} else {
					// You gave us not enough data types to know everything, but
					// some different ones? What are we supposed to do?
					throw new IllegalArgumentException("Size of componentSizes does not match size of components.");
				}
			}
		}
		if (ambiguitiesSizes == null) {
			stateCountsIncludingAmbiguities = new Integer[stateCountsExcludingAmbiguities.length];
			for (int i=0; i<stateCountsExcludingAmbiguities.length; ++i) {
				stateCountsIncludingAmbiguities[i] = stateCountsExcludingAmbiguities[i]
						+ reserveAmbiguities.get().getValue();
			}
		} else {
			stateCountsIncludingAmbiguities = ambiguitiesSizes.getValues();
			for (int i=0; i<stateCountsExcludingAmbiguities.length; ++i) {
				if (stateCountsIncludingAmbiguities[i] < stateCountsExcludingAmbiguities[i]) {
					throw new RuntimeException("Ambiguities cannot reduce the number of different values!");
				}
			}
		}
		System.out.printf("%s derived internal sizes to be %s\n", getID(),
				Arrays.toString(stateCountsIncludingAmbiguities));
	}

	static public int[] compoundState2componentStates(Integer[] components, int compoundState) {
		int[] result = new int[components.length];
		for (int i = components.length - 1; i >= 0; --i) {
			result[i] = compoundState % components[i];
			compoundState /= components[i];
		}
		return result;
	}

	static public int compoundState2componentState(Integer[] components, int compoundState, int component) {
		/*
		 * // Old Implementation: for (int i = components.length - 1; i >
		 * component; --i) { compoundState /= components[i]; } return
		 * compoundState % components[component];
		 */
		return compoundState2componentStates(components, compoundState)[component];
	}

	public int compoundState2componentState(int compoundState, int component) {
		return compoundState2componentState(stateCountsExcludingAmbiguities, compoundState, component);
	}

	static public int componentState2compoundState(Integer[] components, int[] componentStates) {
		int compoundState = 0;
		for (int i = 0; i < components.length; ++i) {
			compoundState *= components[i];
			compoundState += componentStates[i];
		}
		return compoundState;
	}

	public int componentState2compoundState(int[] componentStates) {
		return componentState2compoundState(stateCountsExcludingAmbiguities, componentStates);
	}

	public int getComponentCount() {
		return components.size();
	}

	public Integer[] getStateCounts() {
		return Arrays.copyOf(stateCountsExcludingAmbiguities, stateCountsExcludingAmbiguities.length);
	}

	@Override
	public int getStateCount() {
		return stateCount;
	}

	/**
	 * Convert a sequence represented by a string into a sequence of integers
	 * representing the state for this data type. Ambiguous states should be
	 * represented by integer numbers higher than getStateCount() throws
	 * exception when parsing error occur *
	 */
	@Override
	public List<Integer> string2state(String sSequence) {
		int index = 0;
		int[] subcodes = new int[components.size()];
		for (String code : sSequence.split(";;")) {
			subcodes[index] = components.get(index).string2state(code).get(0);
			++index;
		}
		List<Integer> result = new ArrayList<Integer>(1);
		result.add(componentState2compoundState(subcodes));
		return result;
	}

	/**
	 * Convert an array of states into a sequence represented by a string. This
	 * is the inverse of string2state() throws exception when State cannot be
	 * mapped *
	 */
	@Override
	public String state2string(List<Integer> nStates) {
		throw new RuntimeException("CompoundDataType cannot parse strings");
	}

	@Override
	public String state2string(int[] nStates) {
		throw new RuntimeException("CompoundDataType cannot parse strings");
	}

	/**
	 * returns an array of length getStateCount() containing the (possibly
	 * ambiguous) states that this state represents.
	 */
	@Override
	public boolean[] getStateSet(int iState) {
		boolean[] result = new boolean[] { true };
		for (int i = 0; i < components.size(); ++i) {
			boolean[] byComponent = components.get(i).getStateSet(compoundState2componentState(iState, i));
			boolean[] new_result = new boolean[result.length * byComponent.length];
			for (int j = 0; j < new_result.length; ++j) {
				new_result[j] = result[j / byComponent.length] & byComponent[j % byComponent.length];
			}
			result = new_result;
		}
		return result;
	}

	/**
	 * returns an array with all non-ambiguous states represented by a state.
	 */
	@Override
	public int[] getStatesForCode(int iState) {
		return new int[] { iState };
	}

	@Override
	public boolean isAmbiguousState(int state) {
		for (int i = 0; i < components.size(); ++i) {
			if (components.get(i).isAmbiguousState(compoundState2componentState(state, i))) {
				return true;
			}
		}
		return false;
	}

	/**
	 * true if the class is completely self contained and does not need any
	 * further initialisation. Notable exception: GeneralDataype
	 */
	@Override
	public boolean isStandard() {
		return false;
	}

	@Override
	public String getTypeDescription() {
		String ans = "Compound of data types";
		for (DataType t : components) {
			ans += " <" + t.getTypeDescription() + ">";
		}
		return ans;
	}

	@Override
	public char getChar(int state) {
		return (char) (state + 'A');
	}

	@Override
	public String getCode(int state) {
		String ans = "";
		int i = 0;
		for (DataType t : components) {
			ans += t.getCode(compoundState2componentState(state, i));
			++i;
			if (i < components.size()) {
				ans += ",";
			}
		}
		return ans;
	}

}
