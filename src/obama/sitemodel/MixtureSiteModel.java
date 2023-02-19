package obama.sitemodel;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;

@Description("Site model that is a mixture of various substitution models")
public class MixtureSiteModel extends SiteModelInterface.Base {
	final public Input<RealParameter> weightVectorInput = new Input<>("weights", "mixture weights that determine contribution of each category. "
			+ "Equal weights if not specified.");
	final public Input<RealParameter> rateVectorInput = new Input<>("rates", "mixture rates that specify rate of each category. "
			+ "Equal rates if not specified.");
	final public Input<List<SubstitutionModel>> substModelInput = 
            new Input<>("substModel", "set of substitution models along branches in the beast.tree", new ArrayList<>(), Validate.REQUIRED);

	RealParameter weightVector;
	RealParameter rateVector;
	List<SubstitutionModel> substModels;
	
	@Override
	public void initAndValidate() {
		substModels = substModelInput.get();

		// pick up weights
		weightVector = weightVectorInput.get();
		if (weightVector == null) {
			Double [] weights = new Double[substModels.size()];
			for (int i = 0; i < weights.length; i++) {
				weights[i] = 1.0 / weights.length;
			}
			weightVector = new RealParameter(weights);
		} else if (weightVector.getDimension() != substModels.size()) {
			throw new IllegalArgumentException("dimension of weights (" + weightVector.getDimension() + ") should be the same "
					+ "as number of mixture components (" + substModels.size() + ")");
		}

		// pick up rates
		rateVector = rateVectorInput.get();
		if (rateVector == null) {
			Double [] rates = new Double[substModels.size()];
			for (int i = 0; i < rates.length; i++) {
				rates[i] = 1.0;
			}
			rateVector = new RealParameter(rates);
		} else if (rateVector.getDimension() != substModels.size()) {
			throw new IllegalArgumentException("dimension of rates (" + rateVector.getDimension() + ") should be the same "
					+ "as number of mixture components (" + substModels.size() + ")");
		}
		
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        if (!forceJava) {
        	Log.warning("=============================================================");
        	Log.warning("== WARNING: MixtureSiteModel does not work with BEAGLE     ==");
        	Log.warning("== WARNING: Turn off BEAGLE using the -java flag for BEAST ==");
        	Log.warning("=============================================================");
        }
	}

	
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, int category, double rate,
			double[] matrix) {
    	final double jointBranchRate = getRateForCategory(category, node) * rate;
		substModels.get(category).getTransitionProbabilities(node, startTime, endTime, jointBranchRate, matrix);
	}
	
	@Override
	public boolean integrateAcrossCategories() {
        return true;
	}

	@Override
	public int getCategoryCount() {
		return substModels.size();
	}

	@Override
	public int getCategoryOfSite(int site, Node node) {
        throw new IllegalArgumentException("Integrating across categories");
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		return rateVector.getArrayValue(category);
	}

	@Override
	public double[] getCategoryRates(Node node) {
		return rateVector.getDoubleValues();
	}

	@Override
	public double getProportionForCategory(int category, Node node) {
		return weightVector.getArrayValue(category);
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		return weightVector.getDoubleValues();
	}

}
	