package beast.tools;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import babel.tools.TreeCombiner;
import beast.app.beauti.BeautiDoc;
import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
//import beast.util.TreeParser;

@Description("Grafts nodes into a tree at specified clades with zero branch lengths in order to fossilise clade values "
		+ "in ancestral reconstruction analyses using the TreeSetLikelihood.")
public class TreeFossiliser extends TreeCombiner {
	final public Input<File> cfgFileInput = new Input<>("cfg", "tab separated configuration file containing two columns: "
			+ "column 1: name of taxon\n"
			+ "column 2: a comma separated list of taxa determining MRCA to graft above in source tree (if no constraints have been specified).");
	
	String [] taxonName;
	Set<String> [] subTaxonSets;
	
	// for debugging:
	int found = 0;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		MemoryFriendlyTreeSet srcTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(srcInput.get().getPath(), 0);
		srcTreeSet.reset();
		Tree tree = srcTreeSet.next();
		
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}

		Set<String> taxa = new LinkedHashSet<>();
		for (String taxon : tree.getTaxaNames()) {
			taxa.add(taxon);
		}
		processCfgFile(taxa);

		srcTreeSet.reset();
		int n = taxonName.length;
		while (srcTreeSet.hasNext()) {
			tree = srcTreeSet.next();
			for (int i = 0; i < n; i++) {
				Node src = getMRCA(tree, subTaxonSets[i]);
				Node parent = src.getParent();
				double len = src.getLength();
				// create intermediary node on branch
				double newHeight = src.getHeight();
				
				Node newNode = new Node();
				newNode.setHeight(newHeight);
				newNode.setParent(parent);
				for (int j = 0; j < parent.getChildCount(); j++) {
					if (parent.getChild(j) == src) {
						parent.setChild(j, newNode);
					}
				}
				newNode.addChild(src);
				src.setParent(newNode);
				
				// create new leaf node
				Node leaf = new Node();
				leaf.setID(taxonName[i]);
				Log.warning("Adding " + taxonName[i]);
				leaf.setHeight(newHeight);
				newNode.addChild(leaf);
				leaf.setParent(newNode);
			}
			out.print(tree.getRoot().toNewick());
			out.println(";");
		}
		 
		Log.err("Done!");
		out.close();
	}

	private void processCfgFile(Set<String> taxa) throws IOException {
		
		String cfg = BeautiDoc.load(cfgFileInput.get());
		String [] strs = cfg.split("\n");
		int n = 0;
		for (String str : strs) {
			if (!str.matches("^\\s*$")) {
				n++;
			}
		}
		subTaxonSets = new Set[n];
		taxonName = new String[n];
		int i = 0;
		for (String str : strs) {
			if (!str.matches("^\\s*$")) {
				String [] strs2 = str.split("\t");
				taxonName[i] = strs2[0];
				subTaxonSets[i] = new HashSet<>();
				for (String taxon : strs2[1].split(",")) {
					subTaxonSets[i].add(taxon);					
				}
				i++;
			}
		}
	}




	public static void main(String[] args) throws Exception {
		new Application(new TreeFossiliser(), "Tree Fossiliser", args);
	}

}
