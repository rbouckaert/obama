package obama;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.substitutionmodel.Frequencies;

@Description("PhyloHMM combines hidden markov model (HMM) with tree likelihoods per site")
public class PhyloHMM extends Distribution {
	final public Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("treelikelihood", 
			"treelikelihood, one for each state", 
			new ArrayList<>(), Validate.REQUIRED);
	final public Input<Frequencies> freqsInput = new Input<>("frequencies", "start frequencies for the HMM. Assumed uniform if not specified");
	final public Input<Function> ratesInput = new Input<>("rates", "rates for HMM transition probabilities", Validate.REQUIRED);
	
	public enum hmmAlgorithm {Viterbi, backwardforward};
	final public Input<hmmAlgorithm> hmmAlgorithmInput = new Input<>("HMMAlgorithm", "which HMM algorithm to use "
			+ "", hmmAlgorithm.Viterbi, hmmAlgorithm.values());
	
	final public Input<String> stateLabelsInput = new Input<>("stateLabels", "comma separated list of labels for each of the states in the HMM");
	final public Input<IntegerParameter> stateToOutputMapInput = new Input<>("stateToOutputMap", "map that links HMM states with an output. "
			+ "If not specified, each state is assumed to have a unique output."); 
	
	// threading
    final public Input<Boolean> useThreadsInput = new Input<>("useThreads", "calculated the distributions in parallel using threads (default true)", true);
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

	
	int siteCount, patternCount;
	/** sitePatternIndex[siteCount] maps site index to pattern index **/
	int [] sitePatternIndex;
	int HMMStateCount, HMMOutputCount; 
	/**  HMMpartials[siteCount][HMMStateCount] **/
	double [][] HMMpartials;
	List<TreeLikelihood> likelihoods;
	/** patternLogP[siteCount][HMMOutputCount] **/
	double [][] patternLogP;
	double [][] storedPatternLogP;
	
	Alignment data;
	Function rates;
	
	Frequencies frequencies;
	int [][] maxIndex = null;
	
	String [] stateLabels;
	
	// is the rate matrix for HMM dense or sparse
	boolean isDense = true;
	
	// threading related members
    boolean useThreads;
    int nrOfThreads;
    ExecutorService exec;
    
    IntegerParameter stateToOutputMap;
    int [] map;

	@Override
	public void initAndValidate() {
		likelihoods = likelihoodsInput.get();
		
		if (stateToOutputMapInput.get() == null) {
			HMMStateCount = likelihoods.size();
			HMMOutputCount = HMMStateCount;
			map = new int[HMMStateCount];
			for (int i = 0; i < HMMStateCount; i++) {
				map[i] = i;
			}
		} else {
			stateToOutputMap = stateToOutputMapInput.get();
			HMMStateCount = stateToOutputMap.getDimension();
			Set<Integer> outputs = new LinkedHashSet<>();
			for (Integer i : stateToOutputMap.getValues()) {
				outputs.add(i);
			}
			HMMOutputCount = outputs.size();
			// make sure outputmap contains numbers 0,...,HMMOutputCount-1 (only)
			for (Integer i : outputs) {
				if (i >= HMMStateCount || i < 0) {
					throw new IllegalArgumentException("stateToOutputMap should only all values from 0 to #outputs - 1, (unlike " + i + ")");
				}
			}
			map = new int [HMMStateCount];
			for (int i = 0; i < HMMStateCount; i++) {
				map[i] = stateToOutputMap.getValue(i);
			}
		}
		
		// sanity check
		if (likelihoods.size() != HMMOutputCount) {
			throw new IllegalArgumentException("Number of outputs (" + HMMOutputCount +") does not "
					+ "match number of likelihoods (" + likelihoods.size() +")\n");
		}
		
		if (stateLabelsInput.get() != null) {
			stateLabels = stateLabelsInput.get().split(",");
		} else {
			stateLabels = new String[HMMStateCount];
			for (int i = 0; i < HMMStateCount; i++) {
				stateLabels[i] = "state" + i;
			}
		}
		
		// sanity check
		if (stateLabels.length != HMMStateCount) {
			throw new IllegalArgumentException("Number of labels (" + stateLabels.length +") does not "
					+ "match number of states in HMM (" + HMMStateCount +")\n");
		}
		
		
		rates = ratesInput.get();
		
		// sanity checks
		if (isDense && rates.getDimension() != HMMStateCount * HMMStateCount) {
			throw new IllegalArgumentException("Dimension of rates should be " + (HMMStateCount * HMMStateCount) + " but is " + 
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
		sitePatternIndex = new int[getSiteCount()];
		for (int i = 0; i < getSiteCount(); i++) {
			sitePatternIndex[i] = data.getPatternIndex(i);
		}
		HMMpartials = new double[getSiteCount()][HMMStateCount];
		patternLogP = new double[HMMOutputCount][];
		
		patternCount = data.getPatternCount();
		storedPatternLogP = new double[HMMOutputCount][patternCount];

        useThreads = useThreadsInput.get() && (ProgramStatus.m_nThreads > 1);
		nrOfThreads = useThreads ? ProgramStatus.m_nThreads : 1;
		if (useThreads && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), ProgramStatus.m_nThreads);
		}
		if (useThreads) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}
	
	}
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		
		// collect pattern log likelihoods
		int workAvailable = 0;
        if (useThreads) {
    		for (int i = 0; i < HMMOutputCount; i++) {
    			if (likelihoods.get(i).isDirtyCalculation()) {
	            	workAvailable++;
	            }
	        }
        }
        if (useThreads && workAvailable > 1) {
            calculateUsingThreads(workAvailable);
        } else if (!useThreads || workAvailable > 0) {
    		for (int i = 0; i < HMMOutputCount; i++) {
    			if (likelihoods.get(i).isDirtyCalculation()) {
    				likelihoods.get(i).calculateLogP();
    				patternLogP[i] = likelihoods.get(i).getPatternLogLikelihoods();
    			}
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
	
    CountDownLatch countDown;

    private void calculateUsingThreads(int dirtyDistrs) {
        try {
            countDown = new CountDownLatch(dirtyDistrs);
            // kick off the threads
    		for (int i = 0; i < HMMOutputCount; i++) {
    			if (likelihoods.get(i).isDirtyCalculation()) {
                    exec.execute(new Calculator(likelihoods.get(i)));
    			}
    		}
            countDown.await();
    		for (int i = 0; i < HMMOutputCount; i++) {
    			if (likelihoods.get(i).isDirtyCalculation()) {
    				patternLogP[i] = likelihoods.get(i).getPatternLogLikelihoods();
    			}
    		}
        } catch (RejectedExecutionException | InterruptedException e) {
            useThreads = false;
            Log.err.println("Stop using threads: " + e.getMessage());
        }
    }

    class Calculator implements Runnable {
        Distribution distr;

        Calculator(Distribution distr) {
            this.distr = distr;
        }

        @Override
		public void run() {
            try {
                distr.calculateLogP();
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + distr.getID());
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }

    } // CoreRunnable
    
	void doForward(double[] freqs) {
		double [][] patternP = new double[patternCount][HMMOutputCount];
		double [] patternLogScale = new double[patternCount];

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
		double [] p0 = getHMMpartials()[0];
		double [] P = patternP[sitePatternIndex[0]];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = freqs[i] * P[map[i]];
		}
		
		double logScale = doForwardLoop(patternP);
		
		double totalP = 0;
		double [] p1 = getHMMpartials()[getSiteCount() - 1];
		for (int u = 0; u < HMMStateCount; u++) {
			totalP += p1[u];
		}
		logP = Math.log(totalP) + logScale;
		int [] weights = data.getWeights();
		for (int i = 0; i < patternCount; i++) {
			logP += patternLogScale[i] * weights[i];
		}

	}

	double doForwardLoop(double [][] patternP) {
		double logScale = 0;

		double [] transitionRates = rates.getDoubleValues();
		
		// forward
		for (int i = 1; i < getSiteCount(); i++) {
			
			double [] p0 = getHMMpartials()[i - 1];
			double [] p1 = getHMMpartials()[i];
			
			double [] P = patternP[sitePatternIndex[i]];
			for (int u = 0; u < HMMStateCount; u++) {
				double sum = 0;
				for (int v = 0; v < HMMStateCount; v++) {
					sum += transitionRates[u * HMMStateCount + v] * p0[v];
				}
				p1[u] = sum * P[map[u]];
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
	public void backward() {
		double [] transitionRates = rates.getDoubleValues();
		normalise(getHMMpartials()[getSiteCount() - 1]);
		
		for (int i = getSiteCount()-2; i >= 0; i--) {
			
			double [] p1 = getHMMpartials()[i];
			double [] p0 = getHMMpartials()[i + 1];
			
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
			maxIndex = new int[getSiteCount()][HMMStateCount];
		}
		double [] p0 = getHMMpartials()[0];
		int siteIndex = sitePatternIndex[0];
		for (int i = 0; i < HMMStateCount; i++) {
			p0[i] = Math.log(freqs[i]) + patternLogP[map[i]][siteIndex];
		}
		
		double [] transitionRates = rates.getDoubleValues();
		// to log
		for (int i = 0; i < transitionRates.length; i++) {
			transitionRates[i] = Math.log(transitionRates[i]);
		}
		
		doViterbiLoop(transitionRates);
		
		double [] p1 = getHMMpartials()[getSiteCount() - 1];
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
	
	void doViterbiLoop(double[] transitionRates) {
		double [] p0;
		double [] p1;
		
		// forward
		for (int i = 1; i < getSiteCount(); i++) {
			
			p0 = getHMMpartials()[i - 1];
			p1 = getHMMpartials()[i];
			
			int siteIndex = sitePatternIndex[i];
			for (int u = 0; u < HMMStateCount; u++) {
				double max = Double.NEGATIVE_INFINITY;
				int iMax = -1;
				for (int v = 0; v < HMMStateCount; v++) {
					if (transitionRates[u * HMMStateCount + v] + p0[v] > max) {
						max = transitionRates[u * HMMStateCount + v] + p0[v];
						iMax = v;
					}
				}
				p1[u] = max + patternLogP[map[u]][siteIndex];
				maxIndex[i][u] = iMax;
			}
		}
		
	}


	public void backwardViterbi() {
		// backward
		int [] path = new int[getSiteCount()];
		double max = Double.NEGATIVE_INFINITY;
		int iMax = -1;
		double [] p1 = getHMMpartials()[getSiteCount() - 1];
		for (int v = 0; v < HMMStateCount; v++) {
			if (p1[v] > max) {
				max = p1[v];
				iMax = v;
			}
		}
		path[getSiteCount() - 1] = iMax;
		logP = p1[iMax];
		Arrays.fill(p1,0); p1[iMax] = 1.0;
		
		for (int i = getSiteCount()-2; i >= 0; i--) {
			iMax = maxIndex[i + 1][path[i+1]];
			path[i] = iMax;			
			Arrays.fill(getHMMpartials()[i],0); 
			getHMMpartials()[i][iMax] = 1.0;
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
		
		for (int i = 0; i < HMMOutputCount; i++) {
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
		if (stateToOutputMap != null && stateToOutputMap.somethingIsDirty()) {
			// assume number of states does not change!
			for (int i = 0; i < HMMStateCount; i++) {
				map[i] = stateToOutputMap.getValue(i);
			}
		}
		return super.requiresRecalculation();
	}


	/** return source and target index for the ith rate **/
	public int[] getLink(int i) {
		return new int[]{i % HMMStateCount, i / HMMStateCount};
	}


	public int getSiteCount() {
		return siteCount;
	}


	public double [][] getHMMpartials() {
		return HMMpartials;
	}
	
} // phyloHMM
