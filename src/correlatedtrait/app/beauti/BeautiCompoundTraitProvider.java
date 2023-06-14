package correlatedtrait.app.beauti;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;
import beast.base.parser.PartitionContext;
import beastclassic.app.beauti.TraitDialog;
import beastclassic.evolution.alignment.AlignmentFromTrait;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import correlatedtrait.evolution.alignment.CompoundAlignment4;
import correlatedtrait.evolution.likelihood.CompoundTreeLikelihood;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;



@Description("BEAUti Compound Trait Provider for potentially correlated traits")
public class BeautiCompoundTraitProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc) {
		try {
            List<String> trees = new ArrayList<String>();
            doc.scrubAll(true, false);
            State state = (State) doc.pluginmap.get("state");
            for (StateNode node : state.stateNodeInput.get()) {
                if (node instanceof Tree) { // && ((Tree) node).m_initial.get() != null) {
                    trees.add(BeautiDoc.parsePartition(((Tree) node).getID()));
                }
            }
            TraitDialog dlg = new TraitDialog(doc, trees);
            if (dlg.showDialog("Create new compound trait")) {
            	String tree = dlg.getTree();
            	String name = dlg.getName();
            	PartitionContext context = new PartitionContext(name, name, name, tree);
            	Object o = doc.pluginmap.get(tree);
            	TreeInterface tree_ = (o instanceof TreeInterface) ? (TreeInterface) o : null;
            			
            	if (tree_ == null) {
            		tree_ = (TreeInterface) doc.pluginmap.get("Tree.t:" + tree);
                	if (tree_ == null) {
                		Alert.showMessageDialog(null, "A tree from another partition must be specified");					
                		return null;
                	}
            	}
            	
            	File file = FXUtils.getLoadFile("Tab delimited file containing trait");
        		UserDataType datatype = processTrait(BeautiDoc.load(file), tree_, doc, name);
        		datatype.setID("");
            	if (datatype == null) {
            		// something went wrong importing the data
            		return null;
            	}
            	
            	AlignmentFromTrait alignment = (AlignmentFromTrait) doc.addAlignmentWithSubnet(context, template.get());
            	
            	List<BEASTInterface> list = new ArrayList<BEASTInterface>();
            	list.add(alignment);
            	editAlignment(alignment, doc);
            	return list;
            }
		} catch (Exception e) {
			e.printStackTrace();
		}
        return null;
	}
	
	private UserDataType processTrait(String str, TreeInterface tree, BeautiDoc doc, String id) {
		StringBuilder b = new StringBuilder();
		List<String> taxa = tree.getTaxonset().asStringList();
		Set<String> done = new HashSet<>();
		Set<String> values1 = new HashSet<>();
		Set<String> values2 = new HashSet<>();
		String [] strs = str.split("\n");
		for (String str2 : strs) {
			String [] strs2 = str2.split("\t");
			if (strs2.length >= 3) {
				String taxon = strs2[0];
				if (done.contains(taxon)) {
					Alert.showMessageDialog(null, "Duplicate entry in data file for taxon '" + taxon + "'");
					return null;
				}
				if (!taxa.contains(taxon)) {
					Alert.showMessageDialog(null, "Unrecognised entry in data file '" + taxon + "'");
					return null;
				}
				done.add(taxon);
				b.append(taxon + "=" + strs2[1]+"-" + strs2[2] + ",\n");
				values1.add(strs2[1]);
				values2.add(strs2[2]);
			}
		}
		b.delete(b.length() - 1, b.length());

		if (values1.size() < 2) {
			Alert.showMessageDialog(null, "Data file should contain at least two values for first trait");
			return null;
		}
		if (values2.size() < 2) {
			Alert.showMessageDialog(null, "Data file should contain at least two values for second trait");
			return null;
		}
		StringBuilder codeMap = new StringBuilder();
		boolean i = false;
		for (String v1 : values1) {
			boolean j = false;
			for (String v2 : values2) {
				int v = i ? (j ? 3 : 1) : (j ? 2 : 0); 
				codeMap.append(v1 + "-" + v2 +"=" + v + ",");
				j = true;
			}
			i = true;
		}
		codeMap.delete(codeMap.length() - 1, codeMap.length());
		
		UserDataType dataType = new UserDataType();
		dataType.initByName("states", 4, "codelength", -1, "codeMap", codeMap.toString());
		

		TraitSet traitSet = new TraitSet();
		traitSet.initByName("value", b.toString(), "traitname", "correlated", "taxa", tree.getTaxonset());
		
		CompoundAlignment4 alignment = new CompoundAlignment4();
		alignment.initByName("userDataType", dataType, "traitSet", traitSet);
		alignment.setID(id);
		doc.registerPlugin(dataType);
		doc.registerPlugin(traitSet);
		doc.registerPlugin(alignment);
				
		return dataType;
	}

	@Override
	public int matches(Alignment alignment) {
		for (BEASTInterface output : alignment.getOutputs()) {
			if (output instanceof CompoundTreeLikelihood) {
				return 20;
			}
		}
		return 0;
	}
	
	
	@Override
	public
	void editAlignment(Alignment alignment, BeautiDoc doc) {
		CompoundTraitInputEditor editor = new CompoundTraitInputEditor(doc);
		CompoundTreeLikelihood likelihood = null;
		for (BEASTInterface output : alignment.getOutputs()) {
			if (output instanceof CompoundTreeLikelihood) {
				likelihood = (CompoundTreeLikelihood) output;
				editor.initPanel(likelihood);
		        
				Dialog dlg = new Dialog();
				DialogPane pane = new DialogPane();
				pane.setContent(editor);
				pane.getButtonTypes().add(ButtonType.CLOSE);
				dlg.setDialogPane(pane);
				dlg.setTitle("Discrete trait editor");
		    	pane.setId("DiscreteTraitEditor");
		        dlg.setResizable(true);
		        dlg.showAndWait();
		        
		        try {
			        AlignmentFromTrait traitData = (AlignmentFromTrait) likelihood.dataInput.get();
			        int stateCount = ((UserDataType) traitData.userDataTypeInput.get()).stateCountInput.get();
		        } catch (Exception e) {
					e.printStackTrace();
				}
		        
				return;
			}
		}
	}
	

}
