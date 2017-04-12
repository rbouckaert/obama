package beast.app.beauti;

import javax.swing.Box;

import beast.app.draw.BooleanInputEditor;
import beast.app.draw.EnumInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.sitemodel.OBAMAModelTestSiteModel;
import beast.evolution.substitutionmodel.OBAMAModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.ModelFrequencies;
import beast.evolution.substitutionmodel.NucleotideRevJumpSubstModel;
import beast.evolution.substitutionmodel.SubstitutionModel;

public class OBAMAModelTestInputEditor extends SiteModelInputEditor {
	private static final long serialVersionUID = 1L;

	@Override
    public Class<?> type() {
        return OBAMAModelTestSiteModel.class;
    }
    
	public OBAMAModelTestInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
		OBAMAModelTestSiteModel siteModel = (OBAMAModelTestSiteModel)input.get();
		SubstitutionModel sm = siteModel.substModelInput.get();
		OBAMAModel substModel = (OBAMAModel) sm;
		EnumInputEditor typeEditor = new EnumInputEditor(doc);
		((Box) getComponent(0)).add(typeEditor, 2);
		
//		BooleanInputEditor useExternalFreqs = new BooleanInputEditor(doc);
//		useExternalFreqs.init(substModel.useExternalFreqsInput, substModel, itemNr, bExpandOption, bAddButtons);
//		((Box) getComponent(0)).add(useExternalFreqs, 3);
		validate();
	}
	
}
