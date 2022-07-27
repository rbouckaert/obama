package obama.app.beauti;


import java.util.ArrayList;
import java.util.List;


import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.SiteModelInputEditor;
import javafx.geometry.Insets;
import javafx.scene.Node;
import javafx.scene.Parent;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import obama.sitemodel.OBAMAModelTestSiteModel;
import obama.substitutionmodel.OBAMAModel;
import obama.substitutionmodel.OBAMA_Blosum62;
import obama.substitutionmodel.OBAMA_CpREV;
import obama.substitutionmodel.OBAMA_DCMut;
import obama.substitutionmodel.OBAMA_Dayhoff;
import obama.substitutionmodel.OBAMA_FLU;
import obama.substitutionmodel.OBAMA_HIVb;
import obama.substitutionmodel.OBAMA_HIVw;
import obama.substitutionmodel.OBAMA_JTT;
import obama.substitutionmodel.OBAMA_LG;
import obama.substitutionmodel.OBAMA_MtArt;
import obama.substitutionmodel.OBAMA_MtMam;
import obama.substitutionmodel.OBAMA_MtREV;
import obama.substitutionmodel.OBAMA_RtREV;
import obama.substitutionmodel.OBAMA_VT;
import obama.substitutionmodel.OBAMA_WAG;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;

public class OBAMAModelTestInputEditor extends SiteModelInputEditor {
	static List<EmpiricalSubstitutionModel> availableModels;
	
	@Override
    public Class<?> type() {
        return OBAMAModelTestSiteModel.class;
    }
    
	public OBAMAModelTestInputEditor(BeautiDoc doc) {
		super(doc);
		if (availableModels == null) {		
			availableModels = new ArrayList<>();
			availableModels.add(new OBAMA_Blosum62());
			availableModels.add(new OBAMA_CpREV());
			availableModels.add(new OBAMA_Dayhoff());
			availableModels.add(new OBAMA_DCMut());
			availableModels.add(new OBAMA_FLU());
			availableModels.add(new OBAMA_HIVb());
			availableModels.add(new OBAMA_HIVw());
			availableModels.add(new OBAMA_JTT());
			availableModels.add(new OBAMA_LG());
			availableModels.add(new OBAMA_MtArt());
			availableModels.add(new OBAMA_MtMam());
			availableModels.add(new OBAMA_MtREV());
			availableModels.add(new OBAMA_RtREV());
			availableModels.add(new OBAMA_VT());
			availableModels.add(new OBAMA_WAG());
		}
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
		
		OBAMAModelTestSiteModel siteModel = (OBAMAModelTestSiteModel)input.get();
		SubstitutionModel sm = siteModel.substModelInput.get();
		OBAMAModel substModel = (OBAMAModel) sm;
		List<EmpiricalSubstitutionModel> models = substModel.substModelInput.get();

		VBox box = (VBox)((HBox)pane.getChildren().get(0)).getChildren().get(0);
		for (EmpiricalSubstitutionModel m : availableModels) {
			addCheckBox(m, models, box);
		}
		//validate();
	}

	private void addCheckBox(EmpiricalSubstitutionModel m, List<EmpiricalSubstitutionModel> models, VBox box) {
		String modelName = m.getClass().getSimpleName();
		String modelLabel = modelName;
		if (modelLabel.startsWith("OBAMA_")) {
			modelLabel = modelLabel.substring(6);
		}
	
		CheckBox checkBox = new CheckBox(modelLabel);
		checkBox.setId(modelName);
		
		boolean selected = false;
		for (EmpiricalSubstitutionModel m0 : models) {
			if (m0.getClass() == m.getClass()) {
				selected = true;
				break;
			}
		}
		checkBox.setSelected(selected);
		checkBox.setOnAction(e -> {
			CheckBox b = (CheckBox) e.getSource();
			String label = b.getId();
			setModel(label, b.isSelected());
		});
		checkBox.setPadding(new Insets(5));
		box.getChildren().add(checkBox);
	}

	private void setModel(String label, boolean selected) {
		OBAMAModelTestSiteModel siteModel = (OBAMAModelTestSiteModel) m_input.get();
		SubstitutionModel sm = siteModel.substModelInput.get();
		OBAMAModel substModel = (OBAMAModel) sm;
		List<EmpiricalSubstitutionModel> models = substModel.substModelInput.get();
		if (!selected) {
			for (int i = 0; i< models.size(); i++) {
				EmpiricalSubstitutionModel m = models.get(i);
				if (m.getClass().getSimpleName().equals(label)) {
					models.remove(i);
					return;
				}
			}
		} else {
			// make sure it is not already in the list
			for (int i = 0; i< models.size(); i++) {
				EmpiricalSubstitutionModel m = models.get(i);
				if (m.getClass().getSimpleName().equals(label)) {
					return;
				}
			}
			// add new instance to list
			for (int i = 0; i< availableModels.size(); i++) {
				EmpiricalSubstitutionModel m = availableModels.get(i);
				if (m.getClass().getSimpleName().equals(label)) {
					models.add(m);
					return;
				}
			}
		}
		
	}
	
}
