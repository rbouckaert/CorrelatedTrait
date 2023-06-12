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

package tsa.parameterclone.splitandmerge;

import java.util.HashSet;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;
import beast.base.util.Randomizer;

@Description("Randomly split a group of parameters in two")
@Citation("Huelsenbeck, J.P., Larget, B., Alfaro, M.E., 2004. "
		+ "Bayesian Phylogenetic Model Selection Using Reversible Jump Markov Chain Monte Carlo. "
		+ "Mol Biol Evol 21, 1123-1133. doi:10.1093/molbev/msh123")
public class SplitOperator extends Operator {
	// Inputs that are changed by the operator
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new RealParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", new IntegerParameter(),
			Validate.REQUIRED);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			(IntegerParameter) null);

	Integer maxIndex;

	@Override
	public void initAndValidate() {
		maxIndex = parametersInput.get().getDimension();
		// RealParameter does not implement java.lang.iterable, so we must do
		// the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new RuntimeException(
						"All entries in groupings must be valid indices of parameters");
			}
		}
	}

	/**
	 * Change the parameter and return the log of the Hastings ratio. Split a
	 * class of joined parameters in two.
	 */
	@Override
	public double proposal() {
		// Find the composition of groups, in particular which ones can be
		// split.

		// If only a split can happen, it has probability 1.
		// If splitting and merging can both happen, the split probability
		// is 1/2.
		int nGroups = 0;
		int nGroupsOfSizeAtLeastTwo = 0;
		int[] trueGroupIndices = new int[parametersInput.get().getDimension()];
		int i = 0;
		Integer newIndex = null;
		for (int size : sizesInput.get().getValues()) {
			if (size > 0) {
				++nGroups;
				if (size > 1) {
					trueGroupIndices[nGroupsOfSizeAtLeastTwo] = i;
					++nGroupsOfSizeAtLeastTwo;
				}
			} else {
				if (newIndex == null) {
					newIndex = i;
				}
			}
			++i;
		}
		
		if (newIndex == null) {
			// System.out.printf("Split: Parameter space exhausted\n");
			return Double.NEGATIVE_INFINITY;
		}

		if (nGroupsOfSizeAtLeastTwo < 1) {
			// System.out.printf("Split: No group large enough to split\n");
			return Double.NEGATIVE_INFINITY;
		}

		int rawSplitIndex = Randomizer.nextInt(nGroupsOfSizeAtLeastTwo);
		int splitIndex = trueGroupIndices[rawSplitIndex];

		HashSet<Integer> splitGroup = new HashSet<Integer>();
		i = 0;
		for (int index : groupingsInput.get().getValues()) {
			if (index == splitIndex) {
				splitGroup.add(i);
			}
			++i;
		}

		int newGroupSize;
		int oldGroupSize;
		do {
			// System.out.printf("Split: Generating a valid split...\n");
			// Generate the SPLIT
			newGroupSize = 0;
			oldGroupSize = splitGroup.size();

			// Go through the old list from the end, and either move or keep
			// entries. Note that index 0 is definitely staying in the old
			// group, so we only iterate while index > 0.
			for (int j : splitGroup) {
				if (Randomizer.nextBoolean()) {
					groupingsInput.get().setValue(j, newIndex);
					++newGroupSize;
					--oldGroupSize;
				} else {
					// In case we run into this twice, make sure to reset the
					// split index
					groupingsInput.get().setValue(j, splitIndex);
				}
			}

			// Moving an entry from one group to another means changing the
			// corresponding value in groupings.
			// In theory, there should be a way without rejections.
			// But for correct-before-efficient reasons, just generate any
			// partition, and reject the operator when the partition is trivial.
		} while (newGroupSize == 0 || oldGroupSize == 0);

		double logJacobian = Math.log(newGroupSize + oldGroupSize)
				- Math.log(newGroupSize) - Math.log(oldGroupSize);

		// In order to keep dimensions matched (cf. Green 1995, p. 716), there
		// needs to be a bijection between the pre-image and the image of this
		// operator and its inverse. This is mitigated by a random distortion of
		// the rates, keeping the sum of rates constant.
		// The proposal ration needs to take that into account.
		double rate = parametersInput.get().getValue(splitIndex);
		double mu = Randomizer.uniform(-oldGroupSize * rate, newGroupSize
				* rate);
		// parametersInput.get().log(0, System.out); System.out.println(mu);
		parametersInput.get()
				.setValue(splitIndex, rate + mu / oldGroupSize);
		parametersInput.get().setValue(newIndex, rate - mu / newGroupSize);
		double bijectionDensity = Math
				.log(rate * (oldGroupSize + newGroupSize));

		// Update the group size caches
		sizesInput.get().setValue(newIndex, newGroupSize);
		sizesInput.get().setValue(splitIndex, oldGroupSize);

		// System.out.printf("Split %d into %d\n", splitIndex, newIndex);
		// Now we calculate the Hastings ratio.

		/*
		 * If, after this, only a merge can happen, that merge has probability
		 * 1. But given that we attempt (and reject) the split move in that case
		 * anyway, we leave that factor out. Similarly for the split probability
		 * now.
		 */
		/*
		 * double logMergeProbability = (nGroupsOfSizeAtLeastTwo == 1 &&
		 * newGroupSize == 1 && oldGroupSize == 1) ? 0 : Math.log(0.5); double
		 * logSplitProbability = 0; if (nGroups >= 2) { logSplitProbability =
		 * Math.log(0.5); }
		 */
		
		// The proposal ratio for for a split move is
		// [ P_m(M') 1/(k' nCr 2) ]/[ P_s(M) 1/N(M) 1/(2^(n_i+n_j-1)-1) 1/(q
		// (n_i+n_j)) ]
		// NOTE: The reference states (k nCr 2), but that seems to be a typo. We
		// use k' = k+1 after a split.
		Double p = -Binomial.logChoose(nGroups + 1, 2)
				+ Math.log(nGroupsOfSizeAtLeastTwo)
				+ Math.log(Math.pow(2, newGroupSize + oldGroupSize - 1) - 1)
				+ bijectionDensity + logJacobian;
		// + logMergeProbability
		// - logSplitProbability
		// System.out.printf("Split: %f\n", Math.exp(p));
		return p;
	}
}
