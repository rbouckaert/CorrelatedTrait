package beast.evolution.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.AlignmentFromTrait;

public class AlignmentFromTraitForDPLace extends AlignmentFromTrait {
	final public Input<File> fileInput = new Input<>("file", "file exported from DPlace database containing two features", Validate.REQUIRED);
	final public Input<String> code1Input = new Input<>("code1", "comma separated list of codes to be mapped to 1 for first feature", "1");
	final public Input<String> code2Input = new Input<>("code2", "comma separated list of codes to be mapped to 1 for second feature", "1");

	Map<String, Integer> data = new HashMap<>();

	@Override
	public void initAndValidate() {
		List<String> tabu = new ArrayList<>();
		for (String str : new String [] {"san", "tgu", "wit", "djk", "daf", "afr", "dep", "srm", "egy", "kzh", "tmr", "hat", "lat"}) {
			tabu.add(str);
		}

		Set<String> code1 = new HashSet<>();
		String str = code1Input.get();
		String [] strs = str.split(",");
		for (String str2 : strs) {
			code1.add(str2.trim());
		}
		Set<String> code2 = new HashSet<>();
		str = code2Input.get();
		strs = str.split(",");
		for (String str2 : strs) {
			code2.add(str2.trim());
		}
		
		
		try {
			File file = fileInput.get();
	        BufferedReader fin = new BufferedReader(new FileReader(file));
	        str = null;
	        // eat up header
	        fin.readLine();
	        str = fin.readLine();
	        strs = split(str);
	        int col1 = -1;
	        int col2 = -1;
	        int iso = -1;
	        for (int i = 0; i < strs.length; i++) {
	        	if (strs[i].startsWith("Code")) {
	        		if (col1 >= 0) 
	        			col2 = i;
	        		else 
	        			col1 = i;
	        	}
	        	if (strs[i].matches("ISO.code")) {
	        		iso = i;
	        	}
	        }
	        if (col1 < 0 || col2 < 0 || iso < 0) {
	        	fin.close();
	        	throw new IllegalArgumentException("Could not find columns with Code: and/or ISO code in them. The file does not appear to be exported from DPLace");
	        }
	        while (fin.ready()) {
	            str = fin.readLine();
	            strs = split(str);
	            String s1 = strs[col1];
	            String s2 = strs[col2];
	            String iso0 = strs[iso];
	            if (iso0.length() >= 3) {
	            	if (!s1.equals("NA") && !s2.equals("NA")) {
		            	if (code1.contains(s1)) {
		            		if (code2.contains(s2)) {
		    	            	data.put(iso0, 3);
		            		} else {
		    	            	data.put(iso0, 2);
		            		}
		            	} else {
		            		if (code2.contains(s2)) {
		    	            	data.put(iso0, 1);
		            		} else {
		    	            	data.put(iso0, 0);
		            		}
		            	}
	            	} else {
    	            	data.put(iso0, 15);
	            		//tabu.add(iso0);
	            	}
	            }
	        }
	        fin.close();
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}

		for (String iso0 : tabu) {
        	data.put(iso0, 15);
		}
		
		
		
		
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
        	if (!isTabu(name, tabu)) {
	        	Integer value = data.get(name);
	        	if (value == null) {
	        		throw new IllegalArgumentException("Trait not specified for " + i);
	        	}
	        	List<Integer> values = new ArrayList<>();
	        	values.add(value);
	      		counts.add(values);
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


	static String[] split(String str) {
		List<String> strs = new ArrayList<>();
		int i = 0;
		String str2 = "";
		while (i < str.length()) {
			char c = str.charAt(i);
			if (c == ',' || c== '\t') {
				strs.add(str2);
				str2 = "";
			} else if (c == '"') {
				do {
					str2 += c;
					i++;
					c = str.charAt(i);
				} while (i < str.length() && c != '"');
				str2 += c;
			} else {
				str2 += c;
			}
			i++;
		}
		strs.add(str2);
		return strs.toArray(new String[]{});
	}


	private boolean isTabu(String name, List<String> tabu) {
		for (String s : tabu) {
			if (name.equals(s)) {
				return true;
			}
		}
		return false;
	}

	
}
