package tsa.correlatedcharacters.beauti;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JOptionPane;

import beast.app.beauti.BeautiAlignmentProvider;
import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.util.NexusParser;
import tsa.correlatedcharacters.polycharacter.CompoundAlignment;

public class CompoundAlignmentProvider extends BeautiAlignmentProvider {
	
    /**
     * return new alignment given files
     * @param doc
     * @param files
     * @return
     */
    public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {

        List<BEASTInterface> selectedBEASTObjects = new ArrayList<>();
        for (File file : files) {
            String fileName = file.getName();
			String fileExtension = fileName.substring(fileName.lastIndexOf(".")).toLowerCase();
			Alignment alignment;

			switch (fileExtension) {
				case ".nex":
				case ".nxs":
				case ".nexus":
					NexusParser parser = new NexusParser();
					try {
						parser.parseFile(file);
						if (parser.filteredAlignments.size() > 0) {
							/**
							 * sanity check: make sure the filters do not
							 * overlap
							 **/
							int[] used = new int[parser.m_alignment.getSiteCount()];
							Set<Integer> overlap = new HashSet<>();
							int partitionNr = 1;
							for (Alignment data : parser.filteredAlignments) {
								int[] indices = ((FilteredAlignment) data).indices();
								for (int i : indices) {
									if (used[i] > 0) {
										overlap.add(used[i] * 10000 + partitionNr);
									} else {
										used[i] = partitionNr;
									}
								}
								partitionNr++;
							}
							if (overlap.size() > 0) {
								String overlaps = "<html>Warning: The following partitions overlap:<br/>";
								for (int i : overlap) {
									overlaps += parser.filteredAlignments.get(i / 10000 - 1).getID()
											+ " overlaps with "
											+ parser.filteredAlignments.get(i % 10000 - 1).getID() + "<br/>";
								}
								overlaps += "The first thing you might want to do is delete some of these partitions.</html>";
								JOptionPane.showMessageDialog(null, overlaps);
							}
							/** add alignments **/
							for (Alignment data : parser.filteredAlignments) {
								// sortByTaxonName(data.sequenceInput.get());
								selectedBEASTObjects.add(data);
							}
						} else {
							selectedBEASTObjects.add(parser.m_alignment);
						}
					} catch (Exception ex) {
						ex.printStackTrace();
						JOptionPane.showMessageDialog(null, "Loading of " + fileName + " failed: " + ex.getMessage());
						return null;
					}
					break;

				case ".xml":
					alignment = (Alignment)getXMLData(file);
					selectedBEASTObjects.add(alignment);
					break;

                default:
                    JOptionPane.showMessageDialog(null,
                            "Unsupported sequence file extension.",
                            "Error", JOptionPane.ERROR_MESSAGE);
                    break;
			}
        }
        for (BEASTInterface beastObject : selectedBEASTObjects) {
    		// ensure ID of alignment is unique
			int k = 0;
			String id = beastObject.getID();
			while (doc.pluginmap.containsKey(id)) {
				k++;
				id = beastObject.getID() + k;
			}
			beastObject.setID(id);

			// For each compound alignment, we have a
			// CorrelatedSubstitutionModel, which needs to have its rates
			// defined.
			// So first, we need to introspect the constructed alignment, and
			// see what state counts it contains.

			Integer[] stateCounts = CompoundAlignment.guessSizes((Alignment) beastObject);
			int states = 1;
			int directions = 0;
			for (int count : stateCounts) {
				states *= count;
				directions += count - 1;
			}

			Double[] rates = new Double[states * directions];
			for (int i = 0; i < states * directions; ++i) {
				rates[i] = 0.;
			}
			rates[0] = 1.;
			RealParameter ratesInput = new RealParameter(rates);
			ratesInput.setID("rateValues." + beastObject.getID());
			doc.addPlugin(ratesInput);

			Integer[] groupings = new Integer[states * directions];
			for (int i = 0; i < states * directions; ++i) {
				groupings[i] = 0;
			}
			IntegerParameter groupingsInput = new IntegerParameter(groupings);
			groupingsInput.setID("groupings." + beastObject.getID());
			doc.addPlugin(groupingsInput);

			Integer[] sizes = new Integer[states * directions];
			for (int i = 0; i < states * directions; ++i) {
				sizes[i] = 0;
			}
			sizes[0] = states * directions;
			IntegerParameter sizesInput = new IntegerParameter(sizes);
			sizesInput.setID("sizes." + beastObject.getID());
			doc.addPlugin(sizesInput);

			// sortByTaxonName(((Alignment) beastObject).sequenceInput.get());
			doc.addAlignmentWithSubnet((Alignment) beastObject, template.get());
		}
        return selectedBEASTObjects;
	}
}
