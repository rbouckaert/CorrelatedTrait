package tsa.evolution.tree;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;

@Description("Tree that is selected from a tree set loaded from file")
public class IndexedTreeFromSet extends Tree {
	final public Input<File> treeSetFileInput = new Input<>("treeSetFile", "file containing a tree set in Nexus format", Validate.REQUIRED);
	final public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	
	public Input<IntegerParameter> indexInput = new Input<>("index", "index parameter that points to a single tree in the tree set."
			+ "Only for that tree, the likelihood is calculated", Validate.REQUIRED);

	List<Tree> trees;

	@Override
	public void initAndValidate() {
			// get trees from file
			NexusParser parser = new NexusParser();
			try {
				parser.parseFile(treeSetFileInput.get());
			} catch (IOException e) {
				throw new RuntimeException(e.getMessage());
			}
			trees = parser.trees;
			if (trees == null || trees.size() == 0) {
				throw new IllegalArgumentException("Could not find any trees in nexus file " + treeSetFileInput.get().getName());			
			}
			
			// remove burn-in portion of trees
			int burnin = trees.size() * burninInput.get() / 100;
			if (burnin < 0) {
				burnin = 0;
			}
			if (burnin > trees.size()) {
				throw new IllegalArgumentException("burnin percentage should not exceed 100%");
			}
			Log.warning.println("Removing " + burnin + " trees as burnin");
			for (int i = 0; i < burnin; i++) {
				trees.remove(0);
			}
			
			// sanity check: make sure all taxa in tree are in tree set
			String [] names = getTaxaNames();
			String [] setNames = trees.get(0).getTaxaNames();
			
			for (String id : names) {
				if (indexOf(setNames, id) == -1) {
					throw new IllegalArgumentException("Cannot find taxon " + id + " from tree in list of taxa names "
							+ "in tree set. This may be due to a spelling error, or different case (matching is case sensitive).");
				}
			}

			// remove taxa from tree-set that are not in tree
			Set<String> tabu = new HashSet<>();
			for (String name : setNames) {
				if (indexOf(names, name) == -1) {
					tabu.add(name);
				}
			}
			if (tabu.size() > 0) {
				Log.warning.println("Removing the following taxa from the tree set (because there is no data for them): " + tabu.toString());
				Log.warning.println(names.length + " taxa left.");
				List<Tree> filteredTrees = new ArrayList<>();
				for (Tree tree : trees) {
					int [] internalNodeNr = new int[1];
					internalNodeNr[0] = names.length;
					Tree newTree = new Tree(removeTaxa(tabu, tree.getRoot(), internalNodeNr));
					filteredTrees.add(newTree);
				}
				trees = filteredTrees;
			}
			

			// relabel tips in tree set to match those in tree
			for (Tree tree : trees) {
				relabel(names, tree.getRoot());
			}
		
			this.assignFrom(trees.get(indexInput.get().getValue()));
			super.initAndValidate();

	}
	protected int indexOf(String[] names, String id) {
		for (int i = 0; i < names.length; i++) {
			if (names[i].equals(id)) {
				return i;
			}
 		}
		return -1;
		
	}

	/** change numbers of tip nodes to match those in names **/
	protected void relabel(String [] names, Node node) {
		if (node.isLeaf()) {
			if (node.getID() != null) {
				int nr = indexOf(names, node.getID());
				node.setNr(nr);
			}
		} else {
			for (Node child : node.getChildren()) {
				relabel(names, child);
			}
		}
	}
		
	protected Node removeTaxa(Set<String> taxaToInclude, Node node, int [] internalNodeNr) {
		if (node.isLeaf()) {
			if (taxaToInclude.contains(node.getID())) {
				return null;
			} else {
				return node;
			}
		} else {
			Node left_ = node.getLeft(); 
			Node right_ = node.getRight(); 
			left_ = removeTaxa(taxaToInclude, left_, internalNodeNr);
			right_ = removeTaxa(taxaToInclude, right_, internalNodeNr);
			if (left_ == null && right_ == null) {
				return null;
			}
			if (left_ == null) {
				return right_;
			}
			if (right_ == null) {
				return left_;
			}
			node.removeAllChildren(false);
			node.addChild(left_);
			node.addChild(right_);
			node.setNr(internalNodeNr[0]);
			internalNodeNr[0]++;
			return node;
		}
	}


	@Override
	protected void store() {
		// super.store();
	}
	
	@Override
	public void restore() {
		this.assignFrom(trees.get(indexInput.get().getValue()));
	}
	
	@Override
	protected boolean requiresRecalculation() {
		this.assignFrom(trees.get(indexInput.get().getValue()));
		return super.requiresRecalculation();
	}
}
