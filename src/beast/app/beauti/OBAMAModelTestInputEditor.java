package beast.app.beauti;


import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JPanel;

import beast.app.draw.EnumInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.sitemodel.OBAMAModelTestSiteModel;
import beast.evolution.substitutionmodel.*;

public class OBAMAModelTestInputEditor extends SiteModelInputEditor {
	private static final long serialVersionUID = 1L;
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
		EnumInputEditor typeEditor = new EnumInputEditor(doc);
		((Box) getComponent(0)).add(typeEditor, 2);

		OBAMAModel substModel = (OBAMAModel) sm;
		List<EmpiricalSubstitutionModel> models = substModel.substModelInput.get();

		int k = 3;
		for (EmpiricalSubstitutionModel m : availableModels) {
			addCheckBox(m, models, k++);
		}
//		BooleanInputEditor useExternalFreqs = new BooleanInputEditor(doc);
//		useExternalFreqs.init(substModel.useExternalFreqsInput, substModel, itemNr, bExpandOption, bAddButtons);
//		((Box) getComponent(0)).add(useExternalFreqs, 3);
		validate();
	}

	private void addCheckBox(EmpiricalSubstitutionModel m, List<EmpiricalSubstitutionModel> models, int offset) {
		String modelLabel = m.getClass().getSimpleName();
		JCheckBox checkBox = new JCheckBox(modelLabel);
		
		boolean selected = false;
		for (EmpiricalSubstitutionModel m0 : models) {
			if (m0.getClass() == m.getClass()) {
				selected = true;
				break;
			}
		}
		checkBox.setSelected(selected);
		checkBox.addActionListener(e -> {
			JCheckBox b = (JCheckBox) e.getSource();
			String label = b.getText();
			setModel(label, b.isSelected());
		});
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.X_AXIS));

		panel.add(checkBox);
		panel.add(Box.createHorizontalGlue());
		
		((Box) getComponent(0)).add(panel, offset);
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
