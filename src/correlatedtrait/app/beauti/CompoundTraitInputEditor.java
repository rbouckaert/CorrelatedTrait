package correlatedtrait.app.beauti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.GuessPatternDialog;
import beastfx.app.inputeditor.ListInputEditor;
import beastfx.app.inputeditor.SmallLabel;
import beastfx.app.util.FXUtils;
import correlatedtrait.evolution.likelihood.CompoundTreeLikelihood;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.EventHandler;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.control.TableColumn.CellEditEvent;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.TreeInterface;
import beastclassic.evolution.alignment.AlignmentFromTrait;



public class CompoundTraitInputEditor extends ListInputEditor {
	public CompoundTraitInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> baseType() {
		return CompoundTreeLikelihood.class;
	}

	CompoundTreeLikelihood likelihood;
	TreeInterface tree;
    TraitSet traitSet;
    TextField traitEntry;
    List<String> sTaxa;
    Object[][] tableData;

    public class LocationMap {
		String taxon;
    	String trait;

    	LocationMap(String taxon, String trait) {
    		this.taxon = taxon;
    		this.trait = trait;
    	}
    	
    	public String getTaxon() {
			return taxon;
		}
		public void setTaxon(String taxon) {
			this.taxon = taxon;
		}
		public String getTrait() {
			return trait;
		}
		public void setTrait(String trait) {
			this.trait = trait;
		}
    }
    
    TableView<LocationMap> table;
    ObservableList<LocationMap> taxonMapping;        
    
    public class Data {
    	
    }
    
    UserDataType dataType;

    String m_sPattern = ".*_(..).*";

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr,	ExpandOption bExpandOption, boolean bAddButtons) {
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        this.itemNr = itemNr;
        m_bAddButtons = bAddButtons;
		this.itemNr = itemNr;
		if (itemNr >= 0) {
			likelihood = (CompoundTreeLikelihood) ((ArrayList<?>)input.get()).get(itemNr);
		} else {
			likelihood = (CompoundTreeLikelihood) ((ArrayList<?>)input.get()).get(0);
		}
	}

	public void initPanel(CompoundTreeLikelihood likelihood_) {
		likelihood = likelihood_;
		m_beastObject = likelihood.dataInput.get();
		try {
			m_input = m_beastObject.getInput("traitSet");
		}catch (Exception e) {
			// TODO: handle exception
		}
		
        tree = likelihood.treeInput.get();
        if (tree != null) {
        	Alignment data = likelihood.dataInput.get();
        	if (!(data instanceof AlignmentFromTrait)) {
        		return;
        	}
    		AlignmentFromTrait traitData = (AlignmentFromTrait) data;
            m_input = traitData.traitInput;
            m_beastObject = traitData;
            traitSet = traitData.traitInput.get();
            
            if (traitSet == null) {
                traitSet = new TraitSet();
                String context = BeautiDoc.parsePartition(likelihood.getID());
                traitSet.setID("traitSet." + context);
                try {
                traitSet.initByName("traitname", "discrete",
                        "taxa", tree.getTaxonset(),
                        "value", "");
                m_input.setValue(traitSet, m_beastObject);
                data.initAndValidate();
                } catch (Exception e) {
					// TODO: handle exception
				}
            }
            
            
            dataType = (UserDataType)traitData.userDataTypeInput.get();

            VBox box = FXUtils.newVBox();

            if (traitSet != null) {
                box.getChildren().add(createButtonBox());
                box.getChildren().add(createListBox());

                HBox box2 = FXUtils.newHBox();
                box2.getChildren().add(createTypeList(0));
                box2.getChildren().add(createTypeList(1));
                box.getChildren().add(box2);
            }
            getChildren().add(box);
            validateInput();
            // synchronise with table, useful when taxa have been deleted
            convertTableDataToDataType();
            convertTableDataToTrait();
        }
    } // init


	Set<String> [][] type = new Set[2][2];

    private Node createTypeList(final int index) {
    	String [] strs = dataType.codeMapInput.get().split(",");
    	type[0][0] = new HashSet<>();
    	type[0][1] = new HashSet<>();
    	type[1][0] = new HashSet<>();
    	type[1][1] = new HashSet<>();
    	
    	for (String str: strs) {
    		System.out.println(str);
    		String [] strs2 = str.split("=");
    		if (strs2[1].trim().equals("0")) {
    			String [] strs3 = strs2[0].split("-");
    			type[0][0].add(strs3[0]);
    			type[1][0].add(strs3[1]);
    		} else if (strs2[1].trim().equals("3")) {
    			String [] strs3 = strs2[0].split("-");
    			type[0][1].add(strs3[0]);
    			type[1][1].add(strs3[1]);
    		}
    	}
    	
    	
    	Set<String> all = new HashSet<>();
    	all.addAll(type[index][0]);
    	all.addAll(type[index][1]);
    	String [] all_ = all.toArray(new String[]{});
    	Arrays.sort(all_);

    	VBox box = FXUtils.newVBox();
    	box.getChildren().add(new Label("Categories for trait " + (index + 1) + " "));
    	for (String str : all_) {
    		final CheckBox checkbox = new CheckBox(str);
    		checkbox.setSelected(type[index][1].contains(str));
    		checkbox.setOnAction(e ->
					toggle(index, (CheckBox) e.getSource())
			);
    		box.getChildren().add(checkbox);
    	}
    	
		return box;
	}

	private void toggle(int index, CheckBox value) {
		System.out.println("Toggle " + index + " " + value);
		String str = value.getText();
		if (value.isSelected()) {
			type[index][0].remove(str);
			type[index][1].add(str);
		} else {
			type[index][1].remove(str);
			type[index][0].add(str);
		}
		
		StringBuilder b = new StringBuilder();
		for(String s0 : type[0][0]) {
			for(String s1 : type[1][0]) {
				b.append(s0 + "-" + s1 + "=0,");
			}			
		}
		for(String s0 : type[0][1]) {
			for(String s1 : type[1][0]) {
				b.append(s0 + "-" + s1 + "=1,");
			}			
		}
		for(String s0 : type[0][0]) {
			for(String s1 : type[1][1]) {
				b.append(s0 + "-" + s1 + "=2,");
			}			
		}
		for(String s0 : type[0][1]) {
			for(String s1 : type[1][1]) {
				b.append(s0 + "-" + s1 + "=3,");
			}			
		}
		b.delete(b.length()-1, b.length());
		dataType.codeMapInput.setValue(b.toString(), dataType);
	}

    private TableView createListBox() {
    	try {
    		traitSet.taxaInput.get().initAndValidate();
    		
        	TaxonSet taxa = tree.getTaxonset();
        	taxa.initAndValidate();
        	sTaxa = taxa.asStringList();
    	} catch (Exception e) {
			// TODO: handle exception
            sTaxa = traitSet.taxaInput.get().asStringList();
		}
        String[] columnData = new String[]{"Name", "Trait"};
        taxonMapping = FXCollections.observableArrayList();
        for (String s : sTaxa) {
        	taxonMapping.add(new LocationMap(s, ""));
        }
        // tableData = new Object[sTaxa.size()][2];
        convertTraitToTableData();
        
        
        
        // set up table.
        // special features: background shading of rows
        // custom editor allowing only Date column to be edited.
//      for (Taxon taxonset2 : m_taxonset) {
//      	if (taxonset2 instanceof TaxonSet) {
//		        for (Taxon taxon : ((TaxonSet) taxonset2).taxonsetInput.get()) {
//		            m_lineageset.add(taxon);
//		            m_taxonMap.put(taxon.getID(), taxonset2.getID());
//		            taxonMapping.add(new TaxonMap(taxon.getID(), taxonset2.getID()));
//		        }
//      	}
//      }

      // set up table.
      // special features: background shading of rows
      // custom editor allowing only Date column to be edited.
      table = new TableView<>();        
      table.setPrefWidth(1024);
      table.setEditable(true);
      table.setItems(taxonMapping);

      TableColumn<LocationMap, String> col1 = new TableColumn<>("Taxon");
      col1.setPrefWidth(500);
      col1.setEditable(false);
      col1.setCellValueFactory(
      	    new PropertyValueFactory<LocationMap,String>("Taxon")
      	);
      table.getColumns().add(col1);        

      TableColumn<LocationMap, String> col2 = new TableColumn<>("Trait");
      col2.setPrefWidth(500);
      col2.setEditable(true);
      col2.setCellValueFactory(
      	    new PropertyValueFactory<LocationMap,String>("Trait")
      	);
      col2.setCellFactory(TextFieldTableCell.forTableColumn());
      col2.setOnEditCommit(
              new EventHandler<CellEditEvent<LocationMap, String>>() {
					@Override
					public void handle(CellEditEvent<LocationMap, String> event) {
						String newValue = event.getNewValue();
						LocationMap location = event.getRowValue();
						location.setTrait(newValue);
						convertTableDataToTrait();
  						validateInput();
					}
				}                
          );
      
      table.getColumns().add(col2);        
      
        return table;
    } // createListBox	private Node createListBox() {


    /* synchronise table with data from traitSet Plugin */
    private void convertTraitToTableData() {
        for (int i = 0; i < tableData.length; i++) {
            tableData[i][0] = sTaxa.get(i);
            tableData[i][1] = "";
        }
        String trait = traitSet.traitsInput.get();
        if (trait.trim().length() == 0) {
        	return;
        }
        String[] sTraits = trait.split(",");
        for (String sTrait : sTraits) {
            sTrait = sTrait.replaceAll("\\s+", " ");
            String[] sStrs = sTrait.split("=");
            String value = null;
            if (sStrs.length != 2) {
            	value = "";
                //throw new Exception("could not parse trait: " + sTrait);
            } else {
            	value = sStrs[1].trim();
            }
            String sTaxonID = sStrs[0].trim();
            int iTaxon = sTaxa.indexOf(sTaxonID);
            if (iTaxon < 0) {
            	System.err.println(sTaxonID);
//                throw new Exception("Trait (" + sTaxonID + ") is not a known taxon. Spelling error perhaps?");
            } else {
	            tableData[iTaxon][0] = sTaxonID;
	            tableData[iTaxon][1] = value;
            }
        }

        if (table != null) {
            table.refresh();
        }
    } // convertTraitToTableData

    /**
     * synchronise traitSet Plugin with table data
     */
    private void convertTableDataToTrait() {
        String sTrait = "";
        //Set<String> values = new HashSet<String>(); 
        for (int i = 0; i < tableData.length; i++) {
            sTrait += sTaxa.get(i) + "=" + tableData[i][1];
            if (i < tableData.length - 1) {
                sTrait += ",\n";
            }
        }
        try {
            traitSet.traitsInput.setValue(sTrait, traitSet);
        } catch (Exception e) {
            e.printStackTrace();
        }
        convertTableDataToDataType();
    }

    private void convertTableDataToDataType() {
        List<String> values = new ArrayList<String>(); 
        for (int i = 0; i < tableData.length; i++) {
        	if (tableData[i][1].toString().trim().length() > 0 && !values.contains(tableData[i][1].toString())) {
        		values.add(tableData[i][1].toString());
        	}
        }
        validateInput();
    }

    /**
     * create box with comboboxes for selection units and trait name *
     */
    private HBox createButtonBox() {
        HBox buttonBox = FXUtils.newHBox();

        Label label = new Label("Trait: ");
        //label.setMaximumSize(new Dimension(1024, 20));
        buttonBox.getChildren().add(label);

        traitEntry = new TextField(traitSet.traitNameInput.get());
        traitEntry.setOnKeyReleased(e->{
			try {
				traitSet.traitNameInput.setValue(traitEntry.getText(), traitSet);
			} catch (Exception ex) {
				// TODO: handle exception
			}        	
        });
        // traitEntry.setColumns(12);
        buttonBox.getChildren().add(traitEntry);
        // buttonBox.add(Box.createHorizontalGlue());

        Button guessButton = new Button("Guess");
        guessButton.setId("guess");
        guessButton.setOnAction(e->guess());
        buttonBox.getChildren().add(guessButton);


        m_validateLabel = new SmallLabel("x", "orange");
        m_validateLabel.setVisible(false);
        buttonBox.getChildren().add(m_validateLabel);
        
        return buttonBox;
    } // createButtonBox
    
    
    private void guess() {
        GuessPatternDialog dlg = new GuessPatternDialog(this, m_sPattern);
        //dlg.setName("GuessPatternDialog");
        String sTrait = "";
        switch (dlg.showDialog("Guess traits from taxon names")) {
        case canceled : return;
        case trait: sTrait = dlg.getTrait();
        	break;
        case pattern:
            String sPattern = dlg.getPattern(); 
            try {
                Pattern pattern = Pattern.compile(sPattern);
                for (String sTaxon : sTaxa) {
                    Matcher matcher = pattern.matcher(sTaxon);
                    if (matcher.find()) {
                        String sMatch = matcher.group(1);
                        if (sTrait.length() > 0) {
                            sTrait += ",";
                        }
                        sTrait += sTaxon + "=" + sMatch;
                    }
                    m_sPattern = sPattern;
                }
            } catch (Exception e) {
                return;
            }
            break;
        }
        try {
        	traitSet.traitsInput.setValue(sTrait, traitSet);
        } catch (Exception e) {
			// TODO: handle exception
		}
        convertTraitToTableData();
        convertTableDataToTrait();
        convertTableDataToDataType();
        repaint();
    }
	
	@Override
	public void validateInput() {
		// check all values are specified
		if (tableData == null) {
			return;
		}
        for (int i = 0; i < tableData.length; i++) {
        	if (tableData[i][1].toString().trim().length() == 0) {
        		m_validateLabel.setVisible(true);
        		m_validateLabel.setTooltip(new Tooltip("trait for " + tableData[i][0] + " needs to be specified"));
        		// m_validateLabel.repaint();
        		return;
        	}
        }
		m_validateLabel.setVisible(false);
		super.validateInput();
	}
}
