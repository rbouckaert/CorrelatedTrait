package beast.evolution.likelihood;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import beast.app.BeastMCMC;
import beast.core.BEASTInterface;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.likelihood.ThreadedTreeLikelihood.TreeLikelihoodCaller;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;

public class ThreadedTreeSetLikelihood extends TreeSetLikelihood {
	
	TreeLikelihood [] treelikelihoods;
	int threadCount;
    private ExecutorService pool = null;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();
    double [] treeLogP;
    
    /** private list of likelihoods, to notify framework of TreeLikelihoods being created in initAndValidate() **/
    final private Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("*","",new ArrayList<>());
	
	class MyTreeLikelihood extends TreeLikelihood {
		@Override
		protected boolean requiresRecalculation() {
			super.requiresRecalculation();
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
	}
	

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		treeLogP  = new double[trees.size()];

		TreeLikelihood treelikelihood0 = (TreeLikelihood) treeLikelihoodInput.get();
		Tree tree0 = (Tree) treelikelihood0.treeInput.get();
        
        // set up threading specific stuff
		threadCount = BeastMCMC.m_nThreads;

		treelikelihoods = new TreeLikelihood[threadCount];
    	pool = Executors.newFixedThreadPool(threadCount);
    	int [] boundaries = calcBoundaryPoints(trees.size());
    	for (int i = 0; i < threadCount; i++) {
    		Alignment data = treelikelihood0.dataInput.get();
    		//try {
				treelikelihoods[i] = new MyTreeLikelihood(); // (TreeLikelihood) treelikelihood0.getClass().newInstance();
			//} catch (InstantiationException | IllegalAccessException e) {
			//	e.printStackTrace();
		    //	}
    		treelikelihoods[i].setID(getID() + i);
    		treelikelihoods[i].getOutputs().add(this);
    		likelihoodsInput.get().add(treelikelihoods[i]);

    		treelikelihoods[i].initByName("data", data, 
    				"tree", duplicate((BEASTInterface)tree0, i), 
    				"siteModel", duplicate((BEASTInterface) treelikelihood0.siteModelInput.get(), i), 
    				"branchRateModel", duplicate(treelikelihood0.branchRateModelInput.get(), i), 
    				"useAmbiguities", treelikelihood0.m_useAmbiguities.get(),
                     "scaling" , treelikelihood0.scaling.get() + ""//,
//                     "useJava", new Boolean(treelikelihood0.getInput("useJava")+"")
    				);
    		
    		likelihoodCallers.add(new TreeLikelihoodCaller(treelikelihoods[i], boundaries[i], boundaries[i+1]));
    	}

	}
	
	private int[] calcBoundaryPoints(int n) {
		int [] boundaries = new int[threadCount + 1];
		int range = n / threadCount;
		for (int i = 0; i < threadCount - 1; i++) {
			boundaries[i+1] = range * (i+1);
		}
		boundaries[threadCount] = n;
		return boundaries;
    }


	@Override
	public double calculateLogP() {
        logP = 0;
        
        try {
			pool.invokeAll(likelihoodCallers);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
        
        int n = trees.size();
        // take average over P from treeLogP
        double max = treeLogP[0];
        for (double d : treeLogP) {
        	max = Math.max(d, max);
        }
        double sum = 0;
        for (int i = 0; i < n; i++) {
        	sum += Math.exp(treeLogP[i] - max);
        }
        sum /= n;
        logP = max + Math.log(sum);
        //System.err.println(sum + " " + logP);
        return logP;
	}

	
	@Override
	protected boolean requiresRecalculation() {
		for (TreeLikelihood t : treelikelihoods) {
			t.requiresRecalculation();
		}
		return true;
	}
	
    /** create new instance of src object, connecting all inputs from src object
     * Note if input is a SubstModel, it is duplicated as well.
     * @param src object to be copied
     * @param i index used to extend ID with.
     * @return copy of src object
     */
    private Object duplicate(BEASTInterface src, int i) {
    	if (src == null) { 
    		return null;
    	}
    	BEASTInterface copy;
		try {
			copy = src.getClass().newInstance();
        	copy.setID(src.getID() + "_" + i);
		} catch (InstantiationException | IllegalAccessException e) {
			e.printStackTrace();
			throw new RuntimeException("Programmer error: every object in the model should have a default constructor that is publicly accessible: " + src.getClass().getName());
		}
        for (Input<?> input : src.listInputs()) {
            if (input.get() != null) {
                if (input.get() instanceof List) {
                    // handle lists
                	//((List)copy.getInput(input.getName())).clear();
                    for (Object o : (List<?>) input.get()) {
                        if (o instanceof BEASTInterface) {
                        	// make sure it is not already in the list
                            copy.setInputValue(input.getName(), o);
                        }
                    }
                } else if (input.get() instanceof SubstitutionModel) {
                	// duplicate subst models
                	BEASTInterface substModel = (BEASTInterface) duplicate((BEASTInterface) input.get(), i);
            		copy.setInputValue(input.getName(), substModel);
            	} else {
                    // it is some other value
            		copy.setInputValue(input.getName(), input.get());
            	}
            }
        }
        copy.initAndValidate();
		return copy;
	}
    class TreeLikelihoodCaller implements Callable<Double> {
        private final TreeLikelihood likelihood;
        private final int start, end;

        public TreeLikelihoodCaller(TreeLikelihood likelihood, int start, int end) {
            this.likelihood = likelihood;
            this.start = start;
            this.end = end;
        }

        public Double call() throws Exception {
  		  	try {
  		    	Tree tree0 = (Tree) likelihood.treeInput.get();
  		  		for (int i = start; i < end; i++) {
  		        	tree0.assignFrom(trees.get(i));
  		  			treeLogP[i] = likelihood.calculateLogP();
  		  		}
  		  	} catch (Exception e) {
  		  		System.err.println("Something went wrong in thread (" + start + " to " + end +")");
				e.printStackTrace();
				System.exit(0);
			}
            return 0.0;
        }

    }

}
