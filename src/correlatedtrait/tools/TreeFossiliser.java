package correlatedtrait.tools;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Grafts nodes into a tree at specified clades with zero branch lengths in order to fossilise clade values "
		+ "in ancestral reconstruction analyses using the TreeSetLikelihood.")
public class TreeFossiliser extends beast.base.inference.Runnable {
	final public Input<File> cfgFileInput = new Input<>("cfg", "tab separated configuration file containing two columns: "
			+ "column 1: name of taxon\n"
			+ "column 2: a comma separated list of taxa determining MRCA to graft above in source tree (if no constraints have been specified).");
	final public Input<TreeFile> srcInput = new Input<>("src","source tree (set) file used as skeleton");
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("[[none]]"));

    boolean [] nodesTraversed;
    int nseen;
	
 	
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

	   protected Node getCommonAncestor(Node n1, Node n2) {
	        // assert n1.getTree() == n2.getTree();
	        if( ! nodesTraversed[n1.getNr()] ) {
	            nodesTraversed[n1.getNr()] = true;
	            nseen += 1;
	        }
	        if( ! nodesTraversed[n2.getNr()] ) {
	            nodesTraversed[n2.getNr()] = true;
	            nseen += 1;
	        }
	        while (n1 != n2) {
		        double h1 = n1.getHeight();
		        double h2 = n2.getHeight();
		        if ( h1 < h2 ) {
		            n1 = n1.getParent();
		            if( ! nodesTraversed[n1.getNr()] ) {
		                nodesTraversed[n1.getNr()] = true;
		                nseen += 1;
		            }
		        } else if( h2 < h1 ) {
		            n2 = n2.getParent();
		            if( ! nodesTraversed[n2.getNr()] ) {
		                nodesTraversed[n2.getNr()] = true;
		                nseen += 1;
		            }
		        } else {
		            //zero length branches hell
		            Node n;
		            double b1 = n1.getLength();
		            double b2 = n2.getLength();
		            if( b1 > 0 ) {
		                n = n2;
		            } else { // b1 == 0
		                if( b2 > 0 ) {
		                    n = n1;
		                } else {
		                    // both 0
		                    n = n1;
		                    while( n != null && n != n2 ) {
		                        n = n.getParent();
		                    }
		                    if( n == n2 ) {
		                        // n2 is an ancestor of n1
		                        n = n1;
		                    } else {
		                        // always safe to advance n2
		                        n = n2;
		                    }
		                }
		            }
		            if( n == n1 ) {
	                    n = n1 = n.getParent();
	                } else {
	                    n = n2 = n.getParent();
	                }
		            if( ! nodesTraversed[n.getNr()] ) {
		                nodesTraversed[n.getNr()] = true;
		                nseen += 1;
		            } 
		        }
	        }
	        return n1;
	    }

		protected Node getMRCA(Tree tree, Set<String> taxa) {
			List<Node> leafs = new ArrayList<>();
			for (Node node : tree.getExternalNodes()) {
				if (taxa.contains(node.getID())) {
					leafs.add(node);
				}
			}

	        nodesTraversed = new boolean[tree.getNodeCount()];
	        Node cur = leafs.get(0);

	        for (int k = 1; k < leafs.size(); ++k) {
	            cur = getCommonAncestor(cur, leafs.get(k));
	        }
			return cur;
		}
		



	public static void main(String[] args) throws Exception {
		new Application(new TreeFossiliser(), "Tree Fossiliser", args);
	}

}
