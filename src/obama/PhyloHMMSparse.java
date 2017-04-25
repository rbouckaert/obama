package obama;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import beast.core.Description;
import beast.core.Input;

@Description("Like PhyloHMM, but using sparse transition rate matrix for HMM part")
public class PhyloHMMSparse extends PhyloHMM {
	final public Input<List<Transition>> transitionsInput = new Input<>("transition", "list of state transitions in order "
			+ "corresponding to the rates.", new ArrayList<>());
	
	int [] from;
	int [] to;
	
	public PhyloHMMSparse() {
		isDense = false;
	}
	
	@Override
	public void initAndValidate() {
		List<Transition> transitions = transitionsInput.get();		
		if (transitions.size() != ratesInput.get().getDimension()) {
			throw new IllegalArgumentException("Number of transitions (" + transitions.size()+ ") "
					+ "should be same as dimension of rates (" + ratesInput.get().getDimension() + ")");
		}

		super.initAndValidate();
		
		// set up from/to arrays

		Map<String, Integer> map = new LinkedHashMap<>();
		if (stateLabelsInput.get() != null) {
			for (int i = 0; i < stateLabels.length; i++) {
				map.put(stateLabels[i], i);
			}
		}
		
		for (Transition t : transitionsInput.get()) {
			if (!map.containsKey(t.getFrom())) {
				if (stateLabelsInput.get() != null) {
					throw new IllegalArgumentException("Transition input from (" + t.getFrom() + ") is not in "
							+ "labels. Typo perhaps?");
				}
				map.put(t.getFrom(), map.size());
			}
			if (!map.containsKey(t.getTo())) {
				if (stateLabelsInput.get() != null) {
					throw new IllegalArgumentException("Transition input to (" + t.getTo() + ") is not in "
							+ "labels. Typo perhaps?");
				}
				map.put(t.getTo(), map.size());
			}
		}

		from = new int[transitions.size()];
		to = new int[transitions.size()];
		int k = 0;
		for (Transition t : transitionsInput.get()) {
			from[k] = map.get(t.getFrom());
			to[k] = map.get(t.getTo());
			k++;
		}
		
		
	}
	
	
	/** sparse implementation of forward Viterbi **/
	@Override
	void doViterbi(double[] freqs) {
		if (maxIndex == null) {
			maxIndex = new int[siteCount][HMMStateCount];
		}
		double [] p0 = HMMpartials[0];
		int siteIndex = sitePatternIndex[0];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = Math.log(freqs[i]) + patternLogP[map[i]][siteIndex];
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
			Arrays.fill(p1, Double.NEGATIVE_INFINITY);
			int [] iMax = maxIndex[i];
			Arrays.fill(iMax, -1);
			for (int j = 0; j < transitionRates.length; j++) {
				int u = to[j];
				int v = from[j];
				if (transitionRates[j] + p0[v]> p1[u]) {
					p1[u] = transitionRates[j] + p0[v];
					iMax[u] = v;
				}
			}
			for (int u = 0; u < HMMStateCount; u++) {
				p1[u] += patternLogP[map[u]][siteIndex];
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
	
	@Override
	void doForward(double[] freqs) {
		double [][] patternP = new double[patternCount][HMMOutputCount];
		double [] patternLogScale = new double[patternCount];
		double logScale = 0;

		for (int k = 0; k < patternCount; k++) {
			double max = 0;
			for (int i = 0; i < HMMOutputCount; i++) {
				max = Math.max(max, patternLogP[i][k]);
			}
			patternLogScale[k] = max;
			for (int i = 0; i < HMMOutputCount; i++) {
				patternP[k][i] = Math.exp(patternLogP[i][k] - max);
			}
		}
		
		// initial state
		double [] p0 = HMMpartials[0];
		double [] P = patternP[sitePatternIndex[0]];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = freqs[i] * P[map[i]];
		}
		
		double [] transitionRates = rates.getDoubleValues();
		double [] p1;
		
		// forward
		for (int i = 1; i < siteCount; i++) {
			
			p0 = HMMpartials[i - 1];
			p1 = HMMpartials[i];
			
			P = patternP[sitePatternIndex[i]];
			
			Arrays.fill(p1, 0);
			for (int j = 0; j < transitionRates.length; j++) {
				int u = to[j];
				int v = from[j];
				p1[u] += transitionRates[j] * p0[v];
			}
			for (int u = 0; u < HMMStateCount; u++) {
				p1[u] = p1[u] * P[map[u]];
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
	@Override
	void backward() {
		double [] transitionRates = rates.getDoubleValues();
		normalise(HMMpartials[siteCount - 1]);
		
		double [] sum = new double[HMMStateCount];
		for (int i = siteCount-2; i >= 0; i--) {
			
			double [] p1 = HMMpartials[i];
			double [] p0 = HMMpartials[i + 1];
			
			Arrays.fill(sum, 0);
			for (int j = 0; j < transitionRates.length; j++) {
				int u = to[j];
				int v = from[j];
				sum[u] += transitionRates[j] * p0[v];
			}
			for (int u = 0; u < HMMStateCount; u++) {
				p1[u] *= sum[u];
			}

			normalise(p1);
			
		}
		
	}
	

}