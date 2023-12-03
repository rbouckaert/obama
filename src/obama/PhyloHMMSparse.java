package obama;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;

@Description("Like PhyloHMM, but using sparse transition rate matrix for HMM part")
public class PhyloHMMSparse extends PhyloHMM {
	final public Input<List<Transition>> transitionsInput = new Input<>("transition", "list of state transitions in order "
			+ "corresponding to the rates.", new ArrayList<>());
	
	int [] from;
	int [] to;
		
	public PhyloHMMSparse() {
	}
	
	@Override
	public void initAndValidate() {
		List<Transition> transitions = transitionsInput.get();
		isDense = (transitions.size() == 0);
		if (!isDense) {
			if (transitions.size() != ratesInput.get().getDimension()) {
				throw new IllegalArgumentException("Number of transitions (" + transitions.size()+ ") "
						+ "should be same as dimension of rates (" + ratesInput.get().getDimension() + ")");
			}
		}

		super.initAndValidate();
		
		if (!isDense) {
			setUpTransitions();
		}
	}
	
	// set up from/to arrays based on transitionsInput
	void setUpTransitions() {
		List<Transition> transitions = transitionsInput.get();

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
	void doViterbiLoop(double[] transitionRates) {
		if (isDense) {
			super.doViterbiLoop(transitionRates);
			return;
		}

		// forward
		for (int i = 1; i < getSiteCount(); i++) {
			
			double [] p0 = getHMMpartials()[i - 1];
			double [] p1 = getHMMpartials()[i];
			
			int siteIndex = sitePatternIndex[i];
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
	}
	
	@Override
	double doForwardLoop(double[][] patternP) {
		if (isDense) {
			return super.doForwardLoop(patternP);			
		}
		
		double logScale = 0;
		double [] transitionRates = rates.getDoubleValues();
		
		// forward
		for (int i = 1; i < getSiteCount(); i++) {
			
			double [] p0 = getHMMpartials()[i - 1];
			double [] p1 = getHMMpartials()[i];
			
			double [] P = patternP[sitePatternIndex[i]];
			
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
		return logScale;
	}

	/** calc state distributions in backward sweep **/
	@Override
	public void backward() {
		if (isDense) {
			super.backward();
			return;
		}

		double [] transitionRates = rates.getDoubleValues();
		normalise(getHMMpartials()[getSiteCount() - 1]);
		
		double [] sum = new double[HMMStateCount];
		for (int i = getSiteCount()-2; i >= 0; i--) {
			
			double [] p1 = getHMMpartials()[i];
			double [] p0 = getHMMpartials()[i + 1];
			
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
	

	@Override
	public int[] getLink(int i) {
		if (isDense) {			
			return super.getLink(i);
		}
		return new int[]{from[i], to[i]};
	}
}
