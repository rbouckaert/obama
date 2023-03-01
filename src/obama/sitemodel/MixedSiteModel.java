package obama.sitemodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

@Description("Mixture site model that allows different substitution models for different sites. "
		+ "To be used in combination with MixedTreeLikelihood")
public class MixedSiteModel extends SiteModelInterface.Base {
	final public Input<IntegerParameter> siteModelIndexInput = new Input<>("siteModelIndex", 
			"identifies substitution model for each site. " + 
			"Its length should be number of sites", Validate.REQUIRED);
	final public Input<List<SubstitutionModel>> mixtureComponentInput = 
            new Input<>("component", "pool of substitution models along branches in the beast.tree", new ArrayList<>(), Validate.REQUIRED);
    final public Input<RealParameter> muParameterInput = new Input<>("mutationRate", "mutation rate (defaults to 1.0)");

    protected RealParameter muParameter;
	private IntegerParameter siteModelIndex;
	private List<SubstitutionModel> mixtureComponent;
	private double [][] freqs; 
	
	public MixedSiteModel() {
		substModelInput.setRule(Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		siteModelIndex = siteModelIndexInput.get();
		
		mixtureComponent = initialiseMixtureComponents();
	
		muParameter = muParameterInput.get();
		if (muParameter == null) {
			muParameter = new RealParameter("1.0");
		} 
		if (muParameter.getDimension() != mixtureComponent.size()) {
			muParameter.setDimension(mixtureComponent.size());			
		}
	}

	public void getSiteModelIndex(int [] matrixIndex) {
		int n = siteModelIndex.getDimension();
		if (n != matrixIndex.length) {
			throw new IllegalArgumentException("Expected site model index of length " + matrixIndex.length + " instead of " + n);
		}
		for (int i = 0; i < n; i++) {
			matrixIndex[i] = siteModelIndex.getValue(i);
		}
	}
	
	protected List<SubstitutionModel> initialiseMixtureComponents() {
		return mixtureComponentInput.get();
	}

	public void getTransitionProbabilities(Node node, double startTime, double endTime, int category, double rate,
			double[] matrix) {
    	final double jointBranchRate = /* getRateForCategory(category, node) */ rate * muParameter.getValue(category);
		mixtureComponent.get(category).getTransitionProbabilities(node, startTime, endTime, jointBranchRate, matrix);
	}

	public List<SubstitutionModel> getMixtureComponents() {
		return mixtureComponent;
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
		// if the site model index was updated, we need to recalculate everything 
		if (siteModelIndex.isDirtyCalculation()) {
			Arrays.fill(dirtySiteModels, true);
			return;
		}
		
		for (int i = 0; i < mixtureComponent.size(); i++) {
			dirtySiteModels[i] = 
					((CalculationNode)mixtureComponent.get(i)).isDirtyCalculation()
					|| muParameter.isDirty(i);
		}
	}

	public double[][] getFrequencies() {
		if (freqs == null) {
			freqs = new double[mixtureComponent.size()][];
		}
		for (int i = 0; i < freqs.length; i++) {
			freqs[i] = mixtureComponent.get(i).getFrequencies();
		}
		return freqs;
	}
	
}
