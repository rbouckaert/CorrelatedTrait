/**
 * 
 */
package tsa.correlatedcharacters.polycharacter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.datatype.StandardData;

/**
 * @author gereon
 *
 */
@Description("Class representing a collection of alignment data")
public class CompoundAlignment extends Alignment {
	// Alignment inputs to be ignored:
	// sequence, statecount, dataType, userDataType
	// 'ascertained' stuff
	// TODO: Construct an abstract base class that does not have these.

	// Consider inputs:
	// stripInvariantSitesInput
	// siteWeightsInput
	public Input<Alignment> alignmentInput = new Input<Alignment>("alignment", "The component sites",
			Validate.REQUIRED);
	protected Alignment alignment;

	public CompoundAlignment(Alignment input) {
		super();
		initAndValidate(input);
	}

	public CompoundAlignment() {
		super();
	}

	static public Integer[] guessSizes(Alignment alignment_) {
		Integer[] guessedSizes = new Integer[alignment_.getSiteCount()];
		for (int site = 0; site < alignment_.getSiteCount(); ++site) {
			guessedSizes[site] = 0;
			for (int i : alignment_.getPattern(alignment_.getPatternIndex(site))) {
				if (i >= guessedSizes[site]) {
					guessedSizes[site] = i + 1;
				}
			}
		}
		return guessedSizes;
	}

	private void initAndValidate(Alignment alignment_) {
		alignment = alignment_;

		// Construct or copy the appropriate data type
		CompoundDataType cdt = new CompoundDataType();
		if (userDataTypeInput.get() instanceof CompoundDataType) {
			cdt = (CompoundDataType) userDataTypeInput.get();
		} else if (dataTypeInput.get() == NUCLEOTIDE) {
			// Guess the data type from the data
			Integer[] guessedSizes = guessSizes(alignment);
			List<DataType> components = new ArrayList<DataType>(); 
			for (int i=0; i<guessedSizes.length; ++i) {
				components.add(alignment.getDataType());
			}
			cdt.initByName("components", components, "componentSizesIncludingAmbiguities",
					new IntegerParameter(guessedSizes));
		} else {
			throw new IllegalArgumentException(
					"CompoundAlignment data type is either a CompoundDataType or derived from Alignment and may not be specified otherwise");
		}
		m_dataType = cdt;

		// Given that we take alignments, we don't need to sort, just to check.
		taxaNames = alignment.getTaxaNames();
		// counts, the list of sequences, starts at everything in state 0.
		// stateCounts, the list of stateCount for each sequence, starts being 1
		// everywhere.
		maxStateCount = 1;
		for (int i = 0; i < taxaNames.size(); ++i) {
			ArrayList<Integer> zero = new ArrayList<Integer>();
			zero.add(0);
			counts.add(zero);
			stateCounts.add(cdt.getStateCount());
		}

		for (int taxon_ = 0; taxon_ < alignment.getTaxonCount(); ++taxon_) {
			for (int site_ = 0; site_ < alignment.getSiteCount(); ++site_) {
				counts.get(taxon_).set(0,
						counts.get(taxon_).get(0) * cdt.getStateCounts()[site_] + alignment.getPattern(taxon_, site_));
			}
		}

		maxStateCount = stateCounts.get(0);

		if (maxStateCount != m_dataType.getStateCount()) {
			throw new RuntimeException("Size of data type (" + m_dataType.getStateCount() + ") and of alignments ("
					+ maxStateCount + ") do not match");
		}

		if (siteWeightsInput.get() != null) {
			throw new RuntimeException("No need to specify siteWeights for only one site, did you do something wrong?");
		}

		// grab data from children
		// Sanity check: make sure sequences are of same length

		calcPatterns();
		Log.info.println(toString(false));
	}

	@Override
	public void initAndValidate() {
		initAndValidate(alignmentInput.get());
	}
}
