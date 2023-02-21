package obama.sitemodel;


import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import obama.likelihood.MixtureTreeLikelihood;

@Description("Site model that is a mixture of various substitution models")
public class MixtureSiteModel extends SiteModelInterface.Base {
	final public Input<RealParameter> weightVectorInput = new Input<>("weights", "mixture weights that determine contribution of each mixture component. "
			+ "Equal weights if not specified.");
	final public Input<RealParameter> rateVectorInput = new Input<>("rates", "mixture rates that specify rate of each mixture component. "
			+ "All rates set to 1.0 if not specified.");
	final public Input<List<SubstitutionModel>> mixtureComponentInput = 
            new Input<>("component", "set of substitution models along branches in the beast.tree", new ArrayList<>(), Validate.REQUIRED);
    final public Input<RealParameter> muParameterInput = new Input<>("mutationRate", "mutation rate (defaults to 1.0)");

	protected RealParameter weightVector;
	protected RealParameter rateVector;
	protected RealParameter muParameter;
	protected List<SubstitutionModel> mixtureComponent;
	
	public MixtureSiteModel() {
		substModelInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() {
		mixtureComponent = initialiseMixtureComponents();

		muParameter = muParameterInput.get();
		if (muParameter == null) {
			muParameter = new RealParameter("1.0");
		}
		// pick up weights
		weightVector = weightVectorInput.get();
		if (weightVector == null) {
			Double [] weights = new Double[mixtureComponent.size()];
			for (int i = 0; i < weights.length; i++) {
				weights[i] = 1.0 / weights.length;
			}
			weightVector = new RealParameter(weights);
		} else if (weightVector.getDimension() != mixtureComponent.size()) {
			throw new IllegalArgumentException("dimension of weights (" + weightVector.getDimension() + ") should be the same "
					+ "as number of mixture components (" + mixtureComponent.size() + ")");
		}

		// pick up rates
		rateVector = rateVectorInput.get();
		if (rateVector == null) {
			Double [] rates = new Double[mixtureComponent.size()];
			for (int i = 0; i < rates.length; i++) {
				rates[i] = 1.0;
			}
			rateVector = new RealParameter(rates);
		} else if (rateVector.getDimension() != mixtureComponent.size()) {
			throw new IllegalArgumentException("dimension of rates (" + rateVector.getDimension() + ") should be the same "
					+ "as number of mixture components (" + mixtureComponent.size() + ")");
		}
		
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        if (!forceJava) {
        	Log.warning("=============================================================");
        	Log.warning("== WARNING: MixtureSiteModel does not work with BEAGLE     ==");
        	Log.warning("== WARNING: Turn off BEAGLE using the -java flag for BEAST ==");
        	Log.warning("=============================================================");
        }
        
        
        for (Object o :  getOutputs()) {
        	if (o instanceof TreeLikelihood) {
        		TreeLikelihood tl = (TreeLikelihood) o;
        		tl.implementationInput.setValue(MixtureTreeLikelihood.class.getName(), tl);
        	}
        	if (o instanceof GenericTreeLikelihood && !(o instanceof MixtureTreeLikelihood)) {
        		// ThreadedTreeLikelihood tl = (ThreadedTreeLikelihood) o;
        		// tl.implementationInput.setValue(MixtureTreeLikelihood.class.getName(), tl);
        		throw new IllegalArgumentException("Expected MixtureSiteModel in " + MixtureTreeLikelihood.class.getName() + " not " + o.getClass().getName());
        	}
        }
	}
	
    protected List<SubstitutionModel> initialiseMixtureComponents() {
    	return mixtureComponentInput.get();
    }
	
	public void getTransitionProbabilities(Node node, double startTime, double endTime, int category, double rate,
			double[] matrix) {
    	final double jointBranchRate = getRateForCategory(category, node) * rate * muParameter.getValue();
		mixtureComponent.get(category).getTransitionProbabilities(node, startTime, endTime, jointBranchRate, matrix);
	}
	
	public List<SubstitutionModel> getMixtureComponents() {
		return mixtureComponent;
	}
	
	@Override
	public boolean integrateAcrossCategories() {
        return true;
	}

	@Override
	public int getCategoryCount() {
		return mixtureComponent.size();
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

	@Override	
	public boolean hasImaginaryEigenvectors() {
		return true;
	}
}
	