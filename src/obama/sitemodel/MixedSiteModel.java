package obama.sitemodel;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;

@Description("Mixture site model that allows different substitution models for different sites. "
		+ "To be used in combination with MixedTreeLikelihood")
public class MixedSiteModel extends SiteModelInterface.Base {
	final public Input<IntegerParameter> siteModelIndexInput = new Input<>("siteModelIndex", 
			"identifies substitution model for each site. " + 
			"Its length should be number of sites", Validate.REQUIRED);
	final public Input<List<SubstitutionModel>> mixtureComponentInput = 
            new Input<>("component", "pool of substitution models along branches in the beast.tree", new ArrayList<>(), Validate.REQUIRED);

	IntegerParameter siteModelIndex;
	public void getSiteModelIndex(int [] matrixIndex) {
		int n = siteModelIndex.getDimension();
		if (n != matrixIndex.length) {
			throw new IllegalArgumentException("Expected matrixindex of length " + n + " instead of " + matrixIndex.length);
		}
		for (int i = 0; i < n; i++) {
			matrixIndex[i] = siteModelIndex.getValue(i);
		}
	}
	List<SubstitutionModel> mixtureComponent;
	
	@Override
	public void initAndValidate() {
		siteModelIndex = siteModelIndexInput.get();
		mixtureComponent = mixtureComponentInput.get();
	}

	@Override
	public boolean integrateAcrossCategories() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public int getCategoryCount() {
		return 1;
	}

	@Override
	public int getCategoryOfSite(int site, Node node) {
		return siteModelIndex.getValue(site);
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getCategoryRates(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getProportionForCategory(int category, Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	public void getDirtySiteModels(boolean[] dirtySiteModels) {
		for (int i = 0; i < mixtureComponent.size(); i++) {
			dirtySiteModels[i] = ((CalculationNode)mixtureComponent.get(i)).isDirtyCalculation();
		}
	}
	
}
