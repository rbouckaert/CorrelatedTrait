package tsa;

import java.io.File;
import java.io.IOException;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;

@Description("Likelihood for a set of trees")
public class TreeSetLikelihood extends GenericTreeLikelihood {
	final public Input<File> treeSetFileInput = new Input<>("treeSetFile", "file containing a tree set in Nexus format", Validate.REQUIRED);
	
	List<Tree> trees;
	TreeLikelihood treelikelihood;
	
	
	public TreeSetLikelihood() {
		treeInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() {
		NexusParser parser = new NexusParser();
		try {
			parser.parseFile(treeSetFileInput.get());
		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		}
		trees = parser.trees;
		
		
		Tree tree0 = new Tree();
		tree0.assignFrom(trees.get(0));
		treelikelihood = new TreeLikelihood();
		treelikelihood.initByName("tree", tree0, 
				"data", dataInput.get(),
				"branchratemodel", branchRateModelInput.get(),
				"siteModel", siteModelInput.get()
				);
	
	}
	
	
	@Override
	public double calculateLogP() {
        logP = 0;
    	Tree tree0 = (Tree) treelikelihood.treeInput.get();
        for (Tree tree : trees) {
        	tree0.assignFrom(tree);
        	logP += treelikelihood.calculateLogP();
        }
        return logP;
	}
	
}
