package obama;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.substitutionmodel.Frequencies;

@Description("PhyloHMM combines hidden markov model (HMM) with tree likelihoods per site")
public class PhyloHMM extends Distribution {
	final public Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("treelikelihood", 
			"treelikelihood, one for each state", 
			new ArrayList<>(), Validate.REQUIRED);
	final public Input<Frequencies> freqsInput = new Input<>("frequencies", "start frequencies for the HMM. Assumed uniform if not specified");
	final public Input<RealParameter> ratesInput = new Input<>("rates", "rates for HMM transition probabilities", Validate.REQUIRED);
	
	enum hmmAlgorithm {Viterbi, backwardforward};
	final public Input<hmmAlgorithm> hmmAlgorithmInput = new Input<>("HMMAlgorithm", "which HMM algorithm to use "
			+ "", hmmAlgorithm.Viterbi, hmmAlgorithm.values());
	
	final public Input<String> stateLabelsInput = new Input<>("stateLabels", "comma separated list of labels for each of the states in the HMM");
	
	int siteCount, patternCount;
	/** sitePatternIndex[siteCount] maps site index to pattern index **/
	int [] sitePatternIndex;
	int HMMStateCount; 
	/**  HMMpartials[siteCount][HMMStateCount] **/
	double [][] HMMpartials;
	List<TreeLikelihood> likelihoods;
	/** patternLogP[siteCount][HMMStateCount] **/
	double [][] patternLogP;
	double [][] storedPatternLogP;
	
	Alignment data;
	RealParameter rates;
	
	Frequencies frequencies;
	int [][] maxIndex = null;
	
	String [] stateLabels;
	
	@Override
	public void initAndValidate() {
		likelihoods = likelihoodsInput.get();
		
		HMMStateCount = likelihoods.size();
		if (stateLabelsInput.get() != null) {
			stateLabels = stateLabelsInput.get().split(",");
			if (stateLabels.length != HMMStateCount) {
				throw new IllegalArgumentException("Number of state labels (" + stateLabels.length +") does not "
						+ "match number of likelihoods (" + HMMStateCount +")\n");
			}
		} else {
			stateLabels = new String[HMMStateCount];
			for (int i = 0; i < HMMStateCount; i++) {
				stateLabels[i] = "state" + i;
			}
		}
		rates = ratesInput.get();
		
		// sanity checks
		if (rates.getDimension() != HMMStateCount * HMMStateCount) {
			throw new IllegalArgumentException("Dimension of rates should be " + (HMMStateCount * HMMStateCount) + " but is" + 
					rates.getDimension() + " or perhaps the number of likelihoods should be " + Math.sqrt(rates.getDimension()) +
					" but is " + HMMStateCount);
		}
		
		data = likelihoods.get(0).dataInput.get();
		for (TreeLikelihood likelihood : likelihoods) {
			if (likelihood.dataInput.get() != data) {
				throw new IllegalArgumentException("alignment of likelihoods should be the same, but likelihood " + 
						likelihood.getID() + " has alignment " + likelihood.dataInput.get().getID() + " while " +
						" likelihood " + likelihoods.get(0).getID() + " has alignment " + data.getID());
			}
		}
		
		frequencies = freqsInput.get();
		if (frequencies == null) {
			this.frequencies = new Frequencies();
			RealParameter frequencies = new RealParameter("" + 1.0/HMMStateCount);
			frequencies.setDimension(HMMStateCount);
			this.frequencies.initByName("frequencies", frequencies);
		}
		if (frequencies.getFreqs().length != HMMStateCount) {
			throw new IllegalArgumentException("Frequencies have wrong dimension: is " + frequencies.getFreqs().length +
					" but expected " + HMMStateCount);
		}
		
		// prep site pattern index
		siteCount = data.getSiteCount();
		sitePatternIndex = new int[siteCount];
		for (int i = 0; i < siteCount; i++) {
			sitePatternIndex[i] = data.getPatternIndex(i);
		}
		HMMpartials = new double[siteCount][HMMStateCount];
		patternLogP = new double[HMMStateCount][];
		
		patternCount = data.getPatternCount();
		storedPatternLogP = new double[HMMStateCount][patternCount];
	}
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		
		// collect pattern log likelihoods
		for (int i = 0; i < HMMStateCount; i++) {
			if (likelihoods.get(i).isDirtyCalculation()) {
				likelihoods.get(i).calculateLogP();
				patternLogP[i] = likelihoods.get(i).getPatternLogLikelihoods();
			}
		}
		
		double [] freqs = frequencies.getFreqs();
		
		switch (hmmAlgorithmInput.get()) {
		case Viterbi:
			doViterbi(freqs);
			break;
		case backwardforward:
			doForward(freqs);
			break;
		}
		
		return logP;
	}
	
	
	
	void doForward(double[] freqs) {
		double [][] patternP = new double[patternCount][HMMStateCount];
		double [] patternLogScale = new double[patternCount];
		double logScale = 0;

		for (int k = 0; k < patternCount; k++) {
			double max = 0;
			for (int i = 0; i < HMMStateCount; i++) {
				max = Math.max(max, patternLogP[i][k]);
			}
			patternLogScale[k] = max;
			for (int i = 0; i < HMMStateCount; i++) {
				patternP[k][i] = Math.exp(patternLogP[i][k] - max);
			}
		}
		
		// initial state
		double [] p0 = HMMpartials[0];
		double [] P = patternP[sitePatternIndex[0]];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = freqs[i] * P[i];
		}
		
		double [] transitionRates = rates.getDoubleValues();
		double [] p1;
		
		// forward
		for (int i = 1; i < siteCount; i++) {
			
			p0 = HMMpartials[i - 1];
			p1 = HMMpartials[i];
			
			P = patternP[sitePatternIndex[i]];
			for (int u = 0; u < HMMStateCount; u++) {
				double sum = 0;
				for (int v = 0; v < HMMStateCount; v++) {
					sum += transitionRates[u * HMMStateCount + v] * p0[v];
				}
				p1[u] = sum * P[u];
			}
			// determine scale
			double max = p1[0];
			for (int u = 1; u < HMMStateCount; u++) {
				max = Math.max(max, p1[u]);
			}
			for (int u = 0; u < HMMStateCount; u++) {
				p1[u] /= max;
			}
			logScale += Math.log(max);
		}
		
		
		double totalP = 0;
		p1 = HMMpartials[siteCount - 1];
		for (int u = 0; u < HMMStateCount; u++) {
			totalP += p1[u];
		}
		logP = Math.log(totalP) + logScale;
		int [] weights = data.getWeights();
		for (int i = 0; i < patternCount; i++) {
			logP += patternLogScale[i] * weights[i];
		}
	}

	/** calc state distributions in backward sweep **/
	void backward() {
		double [] transitionRates = rates.getDoubleValues();
		normalise(HMMpartials[siteCount - 1]);
		
		for (int i = siteCount-2; i >= 0; i--) {
			
			double [] p1 = HMMpartials[i];
			double [] p0 = HMMpartials[i + 1];
			
			for (int u = 0; u < HMMStateCount; u++) {
				double sum = 0;
				for (int v = 0; v < HMMStateCount; v++) {
					sum += transitionRates[u * HMMStateCount + v] * p0[v];
				}
				p1[u] *= sum;
			}
			normalise(p1);
			
		}
		
	}

	/** scale to ensure distribution adds to 1 **/
	void normalise(double[] p1) {
		double sum = 0;
		for (double d : p1) {
			sum += d;
		}
		for (int u = 0; u < HMMStateCount; u++) {
			p1[u] /= sum;
		}
	}


	/** dense implementation of forward Viterbi **/
	void doViterbi(double[] freqs) {
		if (maxIndex == null) {
			maxIndex = new int[siteCount][HMMStateCount];
		}
		double [] p0 = HMMpartials[0];
		int siteIndex = sitePatternIndex[0];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = Math.log(freqs[i]) + patternLogP[i][siteIndex];
		}
		
		double [] transitionRates = rates.getDoubleValues();
		// to log
		for (int i = 0; i < transitionRates.length; i++) {
			transitionRates[i] = Math.log(transitionRates[i]);
		}
		double [] p1;
		
		// forward
		for (int i = 1; i < siteCount; i++) {
			
			p0 = HMMpartials[i - 1];
			p1 = HMMpartials[i];
			
			siteIndex = sitePatternIndex[i];
			for (int u = 0; u < HMMStateCount; u++) {
				double max = Double.NEGATIVE_INFINITY;
				int iMax = -1;
				for (int v = 0; v < HMMStateCount; v++) {
					if (transitionRates[u * HMMStateCount + v] + p0[v]> max) {
						max = transitionRates[u * HMMStateCount + v] + p0[v];
						iMax = v;
					}
				}
				p1[u] = max + patternLogP[u][siteIndex];
				maxIndex[i][u] = iMax;
			}
		}
		
		p1 = HMMpartials[siteCount - 1];
		double max = Double.NEGATIVE_INFINITY;
		int iMax = -1;
		for (int v = 0; v < HMMStateCount; v++) {
			if (p1[v] > max) {
				max = p1[v];
				iMax = v;
			}
		}
		logP = p1[iMax];
	}
	
	void backwardViterbi() {
		// backward
		int [] path = new int[siteCount];
		double max = Double.NEGATIVE_INFINITY;
		int iMax = -1;
		double [] p1 = HMMpartials[siteCount - 1];
		for (int v = 0; v < HMMStateCount; v++) {
			if (p1[v] > max) {
				max = p1[v];
				iMax = v;
			}
		}
		path[siteCount - 1] = iMax;
		logP = p1[iMax];
		Arrays.fill(p1,0); p1[iMax] = 1.0;
		
		for (int i = siteCount-2; i >= 0; i--) {
			iMax = maxIndex[i + 1][path[i+1]];
			path[i] = iMax;			
			Arrays.fill(HMMpartials[i],0); 
			HMMpartials[i][iMax] = 1.0;
		}
	}


	@Override
	public List<String> getArguments() {return null;}
	@Override
	public List<String> getConditions() {return null;}
	@Override
	public void sample(State state, Random random) {}
	

	@Override
	public void store() {
		super.store();
		
		for (int i = 0; i < HMMStateCount; i++) {
			System.arraycopy(patternLogP[i], 0, storedPatternLogP[i], 0, patternCount);
		}
	}
	
	@Override
	public void restore() {
		super.restore();
		
		double [][] tmp = patternLogP;
		patternLogP = storedPatternLogP;
		storedPatternLogP = tmp;
	}

	@Override
	protected boolean requiresRecalculation() {
		return super.requiresRecalculation();
	}
	
} // phyloHMM
