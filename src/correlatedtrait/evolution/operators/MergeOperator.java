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

package correlatedtrait.evolution.operators;

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

@Description("Randomly merge two groups of parameters")
@Citation("Huelsenbeck, J.P., Larget, B., Alfaro, M.E., 2004. "
		+ "Bayesian Phylogenetic Model Selection Using Reversible Jump Markov Chain Monte Carlo. "
		+ "Mol Biol Evol 21, 1123-1133. doi:10.1093/molbev/msh123")
public class MergeOperator extends Operator {
	// Inputs that are changed by the operator
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", Validate.REQUIRED);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			(IntegerParameter) null);

	Integer maxIndex;

	@Override
	public void initAndValidate() {
		maxIndex = parametersInput.get().getDimension();
		// Array-like RealParameters do not implement java.lang.iterable, so we
		// must do the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new RuntimeException(
						"All entries in groupings must be valid indices of parameters");
			}
		}
		if (sizesInput.get().getDimension() != parametersInput.get()
				.getDimension()) {
			throw new RuntimeException(
					"sizes must correspond to parameters in dimension");
		}
	}

	/**
	 * Change the parameter and return the log of the Hastings ratio: Merge two
	 * groups of joined parameters, averaging the parameter.
	 */
	@Override
	public double proposal() {

		// Find the composition of groups

		// For each grouping value, the corresponding indices of groupings

		int nGroups = 0;
		int groupsOfSizeAtLeastTwo = 0;
		int[] trueGroupIndices = new int[parametersInput.get().getDimension()];
		int i = 0;
		for (int size : sizesInput.get().getValues()) {
			if (size > 0) {
				trueGroupIndices[nGroups] = i;
				++nGroups;
				if (size > 1) {
					++groupsOfSizeAtLeastTwo;
				}
			}
			++i;
		}

		if (nGroups < 2) {
			// System.out.printf("Merge: No two groups to merge");
			return Double.NEGATIVE_INFINITY;
		}

		int rawMergeIndex = Randomizer.nextInt(nGroups);
		int rawRemoveIndex = Randomizer.nextInt(nGroups - 1);
		// Make sure that the other group is different
		if (rawRemoveIndex >= rawMergeIndex) {
			++rawRemoveIndex;
		}
		int removeIndex = trueGroupIndices[rawRemoveIndex];
		int mergeIndex = trueGroupIndices[rawMergeIndex];

		HashSet<Integer> mergeGroup = new HashSet<Integer>();
		HashSet<Integer> removeGroup = new HashSet<Integer>();
		i = 0;
		for (int index : groupingsInput.get().getValues()) {
			if (index == mergeIndex) {
				mergeGroup.add(i);
			}
			if (index == removeIndex) {
				removeGroup.add(i);
			}
			++i;
		}

		// Generate the MERGE
		Integer mergeGroupSize = mergeGroup.size();
		Integer removeGroupSize = removeGroup.size();

		for (Integer toBeMerged : removeGroup) {
			// groupings[toBeMerged] = mergeIndex
			groupingsInput.get().setValue(toBeMerged, mergeIndex);
		}

		Double logJacobian = Math.log(mergeGroupSize)
				+ Math.log(removeGroupSize)
				- Math.log(mergeGroupSize + removeGroupSize);

		// The merge takes a weighted mean, to conserve the sum of rates.
		double mergedRates = (parametersInput.get().getValue(mergeIndex)
				* mergeGroupSize + parametersInput.get().getValue(
				removeIndex)
				* removeGroupSize)
				/ (mergeGroupSize + removeGroupSize);
		parametersInput.get().setValue(mergeIndex, mergedRates);
		// In order to keep dimensions matched (cf. Green 1995, p. 716), there
		// needs to be a bijection between the pre-image and the image of this
		// operator and its inverse. This is mitigated by a random variable in
		// SplitOperator. The proposal ration needs to take that into account.
		double bijectionDensity = Math.log(mergedRates
				* (mergeGroupSize + removeGroupSize));

		// Update the group size caches
		sizesInput.get().setValue(removeIndex, 0);
		sizesInput.get().setValue(mergeIndex,
				(mergeGroupSize + removeGroupSize));

		// System.out.printf("Merge %d into %d\n", removeIndex, mergeIndex);
		// Now we calculate the Hastings ratio.

		// If we merged two groups of size one, we gain a group of size at
		// least two.
		if (mergeGroupSize == 1 && removeGroupSize == 1) {
			++groupsOfSizeAtLeastTwo;
		}
		// If we merged two groups of size at least two, we lose one in
		// number.
		if (mergeGroupSize >= 2 && removeGroupSize >= 2) {
			--groupsOfSizeAtLeastTwo;
		}

		/*
		 * If, after this, only a split can happen, that split has probability
		 * 1. But given that we attempt the merge move in that case anyway, we
		 * leave that factor out. Similarly for the merge probability now.
		 */
		/*
		 * double logSplitProbability = 0; double logMergeProbability =
		 * groupsOfSizeAtLeastTwo > 0 ? Math.log(0.5) : 0;
		 * 
		 * if (nGroups - 1 >= 2) { logSplitProbability = Math.log(0.5); }
		 */

		// The proposal ratio for for a merge move is
		// [ P_s(M') 1/N(M') 1/(2^(n'_i+n'_j-1)-1) 1/(q' (n'_i+n'_j)) ]/[
		// P_m(M) 1/(k nCr 2) ]

		Double p = -Math.log(groupsOfSizeAtLeastTwo)
				- Math.log(Math.pow(2, mergeGroupSize + removeGroupSize - 1) - 1)
				- bijectionDensity + Binomial.logChoose(nGroups, 2)
				+ logJacobian;
		// + logSplitProbability
		// - logMergeProbability
		// System.out.printf("Merge: %f\n", Math.exp(p));
		return p;

	}
}
