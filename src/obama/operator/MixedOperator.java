package obama.operator;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.BeagleTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import obama.likelihood.MixedTreeLikelihood;
import obama.sitemodel.MixedSiteModel;

@Description("Gibbs operator that proposes new assignment of mixture index")
public class MixedOperator extends Operator {
	final public Input<IntegerParameter> indexInput = new Input<>("index", "index that identifies for each site (or pattern) the component that is used as substitution model", Validate.REQUIRED);
	final public Input<MixedTreeLikelihood> likelihoodInput = new Input<>("likelihood", "mixed tree likelihood used for calculating site likelihoods for each site (or pattern)");
	final public Input<MixedSiteModel> siteModelInput = new Input<>("siteModel", "site model containing mixture components to choose from");

	private IntegerParameter index;
	private MyTreeLikelihood treelikelihood;
	private MyBeagleTreeLikelihood beagleTreeLikelihood;
	private List<SubstitutionModel> components;
	private double [][] siteLogProbs;
	private double [] probs;
	private boolean useSitesNotPatterns;
	private int siteCount;
	private Alignment data;
	
	@Override
	public void initAndValidate() {
		components = siteModelInput.get().getMixtureComponents();
		initTreeLikelihood();
		index = indexInput.get();
		
		int patternCount = likelihoodInput.get().useSitesNotPatternsInput.get() ?
				likelihoodInput.get().dataInput.get().getSiteCount():
				likelihoodInput.get().dataInput.get().getPatternCount();
		int componentCount =  components.size();
		
		siteLogProbs = new double [componentCount][patternCount];
		probs = new double[componentCount];
		
		useSitesNotPatterns = likelihoodInput.get().useSitesNotPatternsInput.get();
		data = likelihoodInput.get().dataInput.get();
		siteCount = data.getSiteCount();
	}

	
	class MyTreeLikelihood extends TreeLikelihood {
		public void setSubstModel(SubstitutionModel substModel, double substitutionRate) {
			this.substitutionModel = substModel;
			this.hasDirt = Tree.IS_FILTHY;
			RealParameter mutationRate = (RealParameter) this.m_siteModel.getInput("mutationRate").get();
			mutationRate.setValue(substitutionRate);
		}
	}
	
	class MyBeagleTreeLikelihood extends BeagleTreeLikelihood {
		public void setSubstModel(SubstitutionModel substModel, double substitutionRate) {
			this.substitutionModel = substModel;
			this.hasDirt = Tree.IS_FILTHY;
			RealParameter mutationRate = (RealParameter) this.m_siteModel.getInput("mutationRate").get();
			mutationRate.setValue(substitutionRate);
		}
	}

	private void initTreeLikelihood() {
		MixedTreeLikelihood mtl = likelihoodInput.get();
		SubstitutionModel dummySubstModel = components.get(0);
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", dummySubstModel, "gammaCategoryCount", 1, "mutationRate", "1.0");
		
		
		MyBeagleTreeLikelihood beagleTreeLikelihood = new MyBeagleTreeLikelihood();
		beagleTreeLikelihood.initByName("data", mtl.dataInput.get(),
                "data", mtl.dataInput.get(), "tree", mtl.treeInput.get(), 
                "siteModel", siteModel,
                "branchRateModel", mtl.branchRateModelInput.get(), 
                "useAmbiguities", mtl.m_useAmbiguities.get(), 
                "useTipLikelihoods", mtl.m_useTipLikelihoods.get(),
                "scaling", mtl.scaling.get().toString());
		if (beagleTreeLikelihood.getBeagle() != null) {
			this.beagleTreeLikelihood = beagleTreeLikelihood;
			return;
		}

		this.beagleTreeLikelihood = null;
		treelikelihood = new MyTreeLikelihood();
		treelikelihood.initByName("data", mtl.dataInput.get(),
                "data", mtl.dataInput.get(), "tree", mtl.treeInput.get(), 
                "siteModel", siteModel,
                "branchRateModel", mtl.branchRateModelInput.get(), 
                "useAmbiguities", mtl.m_useAmbiguities.get(), 
                "useTipLikelihoods", mtl.m_useTipLikelihoods.get(),
                "scaling", mtl.scaling.get().toString());

	}

	@Override
	public double proposal() {
		for (int i = 0; i < siteLogProbs.length; i++) {
			calcLikelihoodForComponent(i);
		}
		
		for (int i = 0; i < siteLogProbs[0].length; i++) {
			int s = sampleIndex(i);
			index.setValue(i, s);
		}		
		
		return Double.POSITIVE_INFINITY;
	}

	private int sampleIndex(int site) {
		for (int i = 0; i < probs.length; i++) {
			probs[i] = siteLogProbs[i][site];
		}
		double max = probs[0];
		for (double d : probs) {
			max = Math.max(max, d);
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] = Math.exp(probs[i] - max);
		}
		double sum = 0;
		for (double d : probs) {
			sum += d;
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		double r = Randomizer.nextDouble();
		int i = 0;
		while (i < probs.length && r > probs[i]) {
			r -= probs[i];
			i++;
		}
		if (i == probs.length) {
			i--;
		}
		return i;
	}

	private void calcLikelihoodForComponent(int i) {
		if (beagleTreeLikelihood != null) {
			beagleTreeLikelihood.setSubstModel(components.get(i), siteModelInput.get().getRateForCategory(i, null));			
			
			beagleTreeLikelihood.calculateLogP();
			
			double [] patternLogLikelihoods = beagleTreeLikelihood.getPatternLogLikelihoods();
			if (useSitesNotPatterns) {
				double [] target = siteLogProbs[i];
				for (int j = 0; j < siteCount; j++) {
					int k = data.getPatternIndex(j);
					target[j] = patternLogLikelihoods[k];
				}
			} else {
				System.arraycopy(patternLogLikelihoods, 0,
					siteLogProbs[i], 0, patternLogLikelihoods.length);
			}
		} else {
			treelikelihood.setSubstModel(components.get(i), siteModelInput.get().getRateForCategory(i, null));
			
			treelikelihood.calculateLogP();
			
			double [] patternLogLikelihoods = treelikelihood.getPatternLogLikelihoods();
			if (useSitesNotPatterns) {
				double [] target = siteLogProbs[i];
				for (int j = 0; j < siteCount; j++) {
					int k = data.getPatternIndex(j);
					target[j] = patternLogLikelihoods[k];
				}
			} else {
				System.arraycopy(patternLogLikelihoods, 0,
					siteLogProbs[i], 0, patternLogLikelihoods.length);
			}
		}
	}

}
