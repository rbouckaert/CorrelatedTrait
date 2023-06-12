package tsa.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.tree.Tree;

@Description("Tree set likelihood that picks only one of the tree at a time")
public class IndexedTreeSetLikelihood extends TreeSetLikelihood {
	public Input<IntegerParameter> indexInput = new Input<>("index", "index parameter that points to a single tree in the tree set."
			+ "Only for that tree, the likelihood is calculated", Validate.REQUIRED);

	
	IntegerParameter index;
	int prevTree;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		index = indexInput.get();
		prevTree = -1;
		index.setLower(0);
		index.setUpper(trees.size());		
	}
	
	@Override
	public double calculateLogP() {
        logP = 0;
    	Tree tree0 = (Tree) treelikelihood.treeInput.get();
    	if (index.getValue() < 0 || index.getValue() >= trees.size()) {
    		logP = Double.NEGATIVE_INFINITY;
    		return logP;
    	}
        tree0.assignFrom(trees.get(index.getValue()));
        tree0.setEverythingDirty(true);
        // flag as dirty to treelikelihood
        // treelikelihood.requiresRecalculation();
        logP += treelikelihood.calculateLogP();
        return logP;
	}
		
}
