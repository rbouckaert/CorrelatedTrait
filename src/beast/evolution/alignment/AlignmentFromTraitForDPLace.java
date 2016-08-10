package beast.evolution.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.evolution.alignment.AlignmentFromTrait;
import beast.evolution.datatype.DataType;
import beast.util.AddOnManager;

public class AlignmentFromTraitForDPLace extends AlignmentFromTrait {
	final static String [] tabu = new String []{"san", "tgu", "wit", "djk", "daf", "afr", "dep", "srm", "egy", "kzh", "tmr", "hat", "lat"};

	
	@Override
	public void initAndValidate() {
    	traitSet = traitInput.get();
    	patternIndex = new int[0];
        counts = new ArrayList<List<Integer>>();
    	if (traitSet == null) { // assume we are in beauti
    		return;
    	}
    	if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
        } else {
        	throw new RuntimeException("expected user data type to be set");
        }

        taxaNames = traitSet.taxaInput.get().taxaNames;
        
        if (traitSet.traitsInput.get() == null || traitSet.traitsInput.get().matches("^\\s*$")) {
        	// prevent initialisation when in beauti
        	patternIndex = new int[1];
            return;
        }
        
        stateCounts = new ArrayList<Integer>();
        for (int i = 0; i < taxaNames.size(); i++) {
        	String name = taxaNames.get(i);
        	if (!isTabu(name)) {
	        	String sValue = traitSet.getStringValue(name);
	        	if (sValue == null) {
	        		throw new IllegalArgumentException("Trait not specified for " + i);
	        	}
	        	List<Integer> iStates = m_dataType.string2state(sValue);
	        	counts.add(iStates);
	        	stateCounts.add(m_dataType.getStateCount());
        	}
        }
        
        List<Taxon> taxa = traitSet.taxaInput.get().taxonsetInput.get();
        for (String s : tabu) {
        	for (int i = 0; i < taxa.size(); i++) {
        		if (taxa.get(i).getID().equals(s)) {
        			taxa.remove(i);
        			break;
        		}
        	}
        }
        traitSet.taxaInput.get().initAndValidate();
        taxaNames = traitSet.taxaInput.get().taxaNames;
        
        
        calcPatterns();
	}


	private boolean isTabu(String name) {
		for (String s : tabu) {
			if (name.equals(s)) {
				return true;
			}
		}
		return false;
	}
}
