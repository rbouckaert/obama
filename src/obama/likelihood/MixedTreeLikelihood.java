/*
* File TreeLikelihood.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/


package obama.likelihood;

import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import obama.sitemodel.MixedSiteModel;

@Description("Tree likelihood that allows site specific site models")
public class MixedTreeLikelihood extends GenericTreeLikelihood {

    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);
    
    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());

    final public Input<Frequencies> rootFrequenciesInput = new Input<>("rootFrequencies", "prior state frequencies at root, optional", Input.Validate.OPTIONAL);

    
    final public Input<Boolean> useSitesNotPatternsInput = new Input<>("useSitesNotPatterns", "Flag to indicate that sites with the same pattern are allowed to use different components. "
    		+ "This is much slower than using patterns, and putting each pattern in the same component", true);

    /**
     * calculation engine *
     */
    protected MixedLikelihoodCore likelihoodCore;
    public MixedLikelihoodCore getLikelihoodCore() {return likelihoodCore;}

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected MixedSiteModel m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}
    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    protected int matrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    protected boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    protected double proportionInvariant = 0;
    public double getProportionInvariant() {return proportionInvariant;}
	public void setProportionInvariant(double proportionInvariant) {this.proportionInvariant = proportionInvariant;
	}
	
	List<Integer> constantPattern = null;
	public List<Integer> getConstantPattern() {return constantPattern;}
    public void setConstantPattern(List<Integer> constantPattern) {this.constantPattern = constantPattern;}
    
    int [] matrixIndex;
    boolean [] dirtySiteModels;

    @Override
    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
		
		if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
			String leaves = "?";
			if (treeInput.get() instanceof Tree) {
				leaves = String.join(", ", ((Tree) treeInput.get()).getTaxaNames());
			}
			throw new IllegalArgumentException(String.format(
					"The number of leaves in the tree (%d) does not match the number of sequences (%d). "
							+ "The tree has leaves [%s], while the data refers to taxa [%s].",
					treeInput.get().getLeafNodeCount(), dataInput.get().getTaxonCount(),
					leaves, String.join(", ", dataInput.get().getTaxaNames())));
		}
        
        
        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (MixedSiteModel) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = getPatternCount();
        likelihoodCore = createLikelihoodCore(stateCount);

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
            if (useSitesNotPatternsInput.get()) {
            	setupAscertainmentSites();
            }
        }
        
        matrixIndex = new int[patterns];
        dirtySiteModels = new boolean[m_siteModel.getMixtureComponents().size()];
    }

    protected MixedLikelihoodCore createLikelihoodCore(int stateCount) {
//		if (stateCount == 4) {
//			return new BeerLikelihoodCore4();
//		} else {
			return new MixedLikelihoodCore(stateCount);
//		}
    }


    /**
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    protected void calcConstantPatternIndices(final int patterns, final int stateCount) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = dataInput.get().getPattern(i);
            final boolean[] isInvariant = new boolean[stateCount];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = dataInput.get().getStateSet(state);
                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousCode(state)) {
                    for (int k = 0; k < stateCount; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < stateCount; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * stateCount + k);
                }
            }
        }
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        
        likelihoodCore.initialize(
                nodeCount,
                getPatternCount(),
                m_siteModel.getMixtureComponents().size(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    private int getPatternCount() {
    	return useSitesNotPatternsInput.get() ? dataInput.get().getSiteCount() : dataInput.get().getPatternCount();
	}
	/**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
	public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
            	int code = useSitesNotPatternsInput.get() ?
                		data.getPattern(taxonIndex, data.getPatternIndex(i))
                		: data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}

	/**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int i = 0; i < patternCount; i++) {
                double[] tipLikelihoods = useSitesNotPatternsInput.get() ?
                		data.getTipLikelihoods(taxonIndex, data.getPatternIndex(i))
                		: data.getTipLikelihoods(taxonIndex, i);
                if (tipLikelihoods != null) {
                	for (int state = 0; state < states; state++) {
                		partials[k++] = tipLikelihoods[state];
                	}
                }
                else {
                	int stateCount = useSitesNotPatternsInput.get() ?
                    		data.getPattern(taxonIndex, data.getPatternIndex(i))
                    		: data.getPattern(taxonIndex, i);
	                boolean[] stateSet = data.getStateSet(stateCount);
	                for (int state = 0; state < states; state++) {
	                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
	                }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    // for testing
    public double[] getRootPartials() {
        return m_fRootPartials.clone();
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
        final TreeInterface tree = treeInput.get();
        
        
		final Node root = tree.getRoot(); // tmp
		for (Node node : root.getAllChildNodesAndSelf()) {
			if (node.isRoot()) continue;
			double length = node.getLength();; 
			if (length < 0) {
				Log.warning("negatuve branch length detected!!! " + length);
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}
		}

        try {
        	m_siteModel.getSiteModelIndex(matrixIndex);
        	m_siteModel.getDirtySiteModels(dirtySiteModels);
        	boolean [] unmask = likelihoodCore.unmask;
        	for (int i = 0; i < matrixIndex.length; i++) {
        		unmask[i] = dirtySiteModels[matrixIndex[i]];
        	}
        	// TODO: remove next line and make recalculation more efficient
        	Arrays.fill(unmask, true);
        	
        	
        	if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
        		calcLogP();
        }
        catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }
    
    
    private List<Integer> m_nIncluded, excludedPatterns;
    private void setupAscertainmentSites() {
    	Alignment data = dataInput.get();
        boolean isAscertained = data.isAscertainedInput.get();

        if (isAscertained) {
            //From AscertainedAlignment
            int from = data.excludefromInput.get();
            int to = data.excludetoInput.get();
            int every = data.excludeeveryInput.get();
            excludedPatterns = new ArrayList<>();
            for (int i = from; i < to; i += every) {
                excludedPatterns.add(i);
            }
            
    		from = data.m_includefrom.get();
    		to = data.m_includeto.get();
    		every = data.m_includeevery.get();
    		m_nIncluded = new ArrayList<>();
    		for (int i = from; i < to; i += every) {
    			m_nIncluded.add(i);
    		}
        }
    } // initAndValidate

    
    private double getAscertainmentCorrection(double[] patternLogProbs) {
        double excludeProb = 0, includeProb = 0, returnProb = 1.0;

        for (int i = 0; i < m_nIncluded.size(); i++) {
        	includeProb += Math.exp(patternLogProbs[m_nIncluded.get(i)]);
        }

        for (int i : excludedPatterns) {
            excludeProb += Math.exp(patternLogProbs[i]);
        }

        if (includeProb == 0.0) {
            returnProb -= excludeProb;
        } else if (excludeProb == 0.0) {
            returnProb = includeProb;
        } else {
            returnProb = 1.0 + includeProb - excludeProb;
        }
        return Math.log(returnProb);
    } // getAscertainmentCorrection

    protected void calcLogP() {
        logP = 0.0;
        if (useSitesNotPatternsInput.get()) {
	        if (useAscertainedSitePatterns) {
	            final double ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
	            for (int i = 0; i < dataInput.get().getSiteCount(); i++) {
	            	int j = dataInput.get().getPatternIndex(i);
	            	if (dataInput.get().getPatternWeight(j) != 0) {
	            		logP += (patternLogLikelihoods[i] - ascertainmentCorrection);
	            	}
	            }
	        } else {
	            for (int i = 0; i < dataInput.get().getSiteCount(); i++) {
	            	int j = dataInput.get().getPatternIndex(i);
	            	if (dataInput.get().getPatternWeight(j) != 0) {
	            		logP += patternLogLikelihoods[i];
	            	}
	            }
	        }        	
        } else {
	        if (useAscertainedSitePatterns) {
	            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
	            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
	                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
	            }
	        } else {
	            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
	                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
	            }
	        }
        }
    }

    /* Assumes there IS a branch rate model as opposed to traverse() */
    protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            if (branchTime != m_branchLengths[nodeIndex]) {
            	Arrays.fill(likelihoodCore.unmask, true);
            }
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            List<SubstitutionModel> components = m_siteModel.getMixtureComponents();
            for (int i = 0; i < components.size(); i++) {
            	// TODO: check component is dirty + branch rate is dirty, otherwise do no need to update this particular matrix 
                m_siteModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), i, branchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                //if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex, matrixIndex);
                //} else {
                //    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                //}

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

//                    final double[] proportions = m_siteModel.getCategoryProportions(node);
//                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

//                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        proportionInvariant = m_siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
//                    }

                    double[][] rootFrequencies = m_siteModel.getFrequencies();
//                    if (rootFrequenciesInput.get() != null) {
//                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
//                    }
                    likelihoodCore.calculateLogLikelihoods(node.getNr(), rootFrequencies, patternLogLikelihoods, matrixIndex);
                }

            }
        }
        return update;
    } // traverseWithBRM

    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
	public double [] getPatternLogLikelihoods() {
		return patternLogLikelihoods.clone();
	} // getPatternLogLikelihoods

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (rootFrequenciesInput.get() != null && rootFrequenciesInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
	public List<String> getArguments() {
        return Collections.singletonList(dataInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
	public List<String> getConditions() {
        return m_siteModel.getConditions();
    }

} // class TreeLikelihood
