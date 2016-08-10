package beast.evolution.likelihood;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.State;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;

@Description("Likelihood for a set of trees")
public class TreeSetLikelihood extends Distribution {
	final public Input<File> treeSetFileInput = new Input<>("treeSetFile", "file containing a tree set in Nexus format", Validate.REQUIRED);
	final public Input<GenericTreeLikelihood> treeLikelihoodInput = new Input<>("treeLikelihood", "tree likelihood used to evaluate each of the trees", Validate.REQUIRED);
	final public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	
	List<Tree> trees;
	TreeLikelihood treelikelihood;
	
	
	public TreeSetLikelihood() {
	}
	
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
		treelikelihood = (TreeLikelihood) treeLikelihoodInput.get();
		Tree tree0 = (Tree) treelikelihood.treeInput.get();
		String [] names = tree0.getTaxaNames();
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
	
	}
	private Node removeTaxa(Set<String> taxaToInclude, Node node, int [] internalNodeNr) {
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


	/** change numbers of tip nodes to match those in names **/
	private void relabel(String [] names, Node node) {
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

	private int indexOf(String[] names, String id) {
		for (int i = 0; i < names.length; i++) {
			if (names[i].equals(id)) {
				return i;
			}
 		}
		return -1;
		
	}

	@Override
	public double calculateLogP() {
        logP = 0;
        int n = trees.size();
        double [] treeLogP  = new double[n];
    	Tree tree0 = (Tree) treelikelihood.treeInput.get();
    	int i = 0;
        for (Tree tree : trees) {
        	tree0.assignFrom(tree);
        	treeLogP[i++] += treelikelihood.calculateLogP();
        }
        
        // take average over P from treeLogP
        double max = treeLogP[0];
        for (double d : treeLogP) {
        	max = Math.max(d, max);
        }
        double sum = 0;
        for (i = 0; i < n; i++) {
        	sum += Math.exp(treeLogP[i] - max);
        }
        sum /= n;
        logP = max + Math.log(sum);
        //System.err.println(sum + " " + logP);
        return logP;
	}

	
	@Override
	protected boolean requiresRecalculation() {
		treelikelihood.requiresRecalculation();
		return true;
	}
	
	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}
	
}
