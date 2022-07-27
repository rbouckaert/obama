package beast.evolution.substitutionmodel;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

@Description("Substitution model for amino acid based on nucleotide models for individual codon positions")
public class ParameterisedAminoAcidModel extends GeneralSubstitutionModel {
	final public Input<GeneralSubstitutionModel> substModel1Input = new Input<>("model1", "nucleotide substitution model for "
			+ "codon position 1", Validate.REQUIRED);
	final public Input<GeneralSubstitutionModel> substModel2Input = new Input<>("model2", "nucleotide substitution model for "
			+ "codon position 2", Validate.REQUIRED);
	final public Input<GeneralSubstitutionModel> substModel3Input = new Input<>("model3", "nucleotide substitution model for "
			+ "codon position 3", Validate.REQUIRED);
	final public Input<RealParameter> substRate1Input = new Input<>("rate1", "substitution rate for codon position 1", new RealParameter("1.0"));
	final public Input<RealParameter> substRate2Input = new Input<>("rate2", "substitution rate for codon position 2", new RealParameter("1.0"));
	final public Input<RealParameter> substRate3Input = new Input<>("rate3", "substitution rate for codon position 3", new RealParameter("1.0"));
	
    public ParameterisedAminoAcidModel() {
        frequenciesInput.setRule(Validate.OPTIONAL);
        ratesInput.setRule(Validate.OPTIONAL);
    }


    RealParameter substRate1, substRate2, substRate3;
    GeneralSubstitutionModel model1, model2, model3;
    List<int []> ratemap1;
    List<int []> ratemap2;
    List<int []> ratemap3;

	final String aacode = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
	final String encoding = new Aminoacid().getCodeMap(); // "ACDEFGHIKLMNPQRSTVWY";
	int [] index;
	
    @Override
    public void initAndValidate() {
    	substRate1 = substRate1Input.get();
    	substRate2 = substRate2Input.get();
    	substRate3 = substRate3Input.get();
    	model1 = substModel1Input.get();
    	model2 = substModel2Input.get();
    	model3 = substModel3Input.get();

    	updateMatrix = true;
        nrOfStates = 20;

        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e.getMessage());
		}

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)];
        
        setupRateMaps();
        setupIndex();
        
        
        if (frequenciesInput.get() != null) {
        	frequencies = frequenciesInput.get();
        } else {
        	
        	
	        String valuesString = "";
	        for (int i = 0; i < nrOfStates; i++) {
	            valuesString += (1.0/nrOfStates) + " ";
	        }
	        RealParameter freqsRParam = new RealParameter(valuesString);
	        frequencies = new Frequencies();
	        frequencies.initByName("frequencies", freqsRParam);
        }
    }




    private void setupIndex() {
    	index = new int[64];
    	for (int i = 0; i < 64; i++) {
    		int aa1 = aacode.charAt(i);
    		int s1 = encoding.indexOf(aa1);
    		index[i] = s1;
    	}
		
	}




	private void setupRateMaps() {
    	ratemap1 = new ArrayList<>();
    	ratemap2 = new ArrayList<>();
    	ratemap3 = new ArrayList<>();
		

    	char aa1, aa2;
		int cs1, cs2, s1, s2;

		for (int i1 = 0; i1 < 4; i1++) {
			for (int j1 = 0; j1 < 4; j1++) {
				for (int k1 = 0; k1 < 4; k1++) {
					cs1 = 16 * k1 + 4 * j1 + i1;
					aa1 = aacode.charAt(cs1);
					s1 = encoding.indexOf(aa1);

					if (s1 >= 0 && s1 < nrOfStates) {
						for (int i2 = 0; i2 < 4; i2++) {
							for (int j2 = 0; j2 < 4; j2++) {
								for (int k2 = 0; k2 < 4; k2++) {
									cs2 = 16 * k2 + 4 * j2 + i2;
									aa2 = aacode.charAt(cs2);
									s2 = encoding.indexOf(aa2);
									
									if (s2 >= 0  && s1 < nrOfStates && s1 != s2) {
										int target = s1 * 19 + s2 + (s2 > s1 ? -1 : 0);
										// only consider one codon position changes
										if (i1 != i2 && j1 == j2 && k1 == k2) {
											ratemap1.add(new int[]{target, i1, i2});
										}
										if (i1 == i2 && j1 != j2 && k1 == k2) {
											ratemap2.add(new int[]{target, j1, j2});
										}
										if (i1 == i2 && j1 == j2 && k1 != k2) {
											ratemap3.add(new int[]{target, k1, k2});
										}
									}
								}
							}
						}
					}
				}
			}
		}
		

		for (int [] i : ratemap1) {
			for (int [] j : ratemap2) {
				if (i[0] == j[0]) {
					System.out.println(Arrays.toString(i) + " " + Arrays.toString(j));
				}
			}
			for (int [] j : ratemap3) {
				if (i[0] == j[0]) {
					System.out.println(Arrays.toString(i) + " " + Arrays.toString(j));
				}
			}
		}
		for (int [] i : ratemap2) {
			for (int [] j : ratemap3) {
				if (i[0] == j[0]) {
					System.out.println(Arrays.toString(i) + " " + Arrays.toString(j));
				}
			}
		}
	}
        
    @Override
    public void setupRelativeRates() {
    	Arrays.fill(relativeRates, 0);

    	model1.setupRelativeRates();
		model1.setupRateMatrix();
    	model2.setupRelativeRates();
		model2.setupRateMatrix();
    	model3.setupRelativeRates();
		model3.setupRateMatrix();
		double [][] rates1 = model1.getRateMatrix();
		double [][] rates2 = model2.getRateMatrix();
		double [][] rates3 = model3.getRateMatrix();
		double rate1 = substRate1.getValue();
		double rate2 = substRate2.getValue();
		double rate3 = substRate3.getValue();
    	
    	char aa1, aa2;
		int cs1, cs2, s1, s2;

		for (int i1 = 0; i1 < 4; i1++) {
			for (int j1 = 0; j1 < 4; j1++) {
				for (int k1 = 0; k1 < 4; k1++) {
					cs1 = 16 * k1 + 4 * j1 + i1;
					aa1 = aacode.charAt(cs1);
					s1 = encoding.indexOf(aa1);

					if (s1 >= 0 && s1 < nrOfStates) {
						for (int i2 = 0; i2 < 4; i2++) {
							for (int j2 = 0; j2 < 4; j2++) {
								for (int k2 = 0; k2 < 4; k2++) {
									cs2 = 16 * k2 + 4 * j2 + i2;
									aa2 = aacode.charAt(cs2);
									s2 = encoding.indexOf(aa2);
									
									if (s2 >= 0  && s1 < nrOfStates && s1 != s2) {
										relativeRates[s1 * 19 + s2 + (s2 > s1 ? -1 : 0)] +=
												(i1 != i2 ? rates1[i1][i2] * rate1: 1.0) * 
												(j1 != j2 ? rates2[j1][j2] * rate2: 1.0) * 
												(k1 != k2 ? rates3[k1][k2] * rate3: 1.0);
									}
								}
							}
						}
					}
				}
			}
		}

//    	processRates(model1, ratemap1, substRate1);
//    	processRates(model2, ratemap2, substRate2);
//    	processRates(model3, ratemap3, substRate3);
    	
    	// is symmetric?
    	for (int i = 0; i < nrOfStates; i++) {
    		for (int j = i + 1; j < nrOfStates; j++) {
    			if (Math.abs(relativeRates[i*19+j-1] - relativeRates[j*19+i]) > 1e-13) {
    				int h = 3;
    				h++;
    			}
    		}
    	}
    }
    
    private void processRates(GeneralSubstitutionModel model1, List<int[]> ratemap1, RealParameter substRate1) {
    	//if (model1.updateMatrix) {
    		model1.setupRelativeRates();
    		model1.setupRateMatrix();
    	//	model1.updateMatrix = false;
    	//}
    	double [][] rates1 = model1.getRateMatrix();
    	double [] relrates1 = model1.getRelativeRates();
    	double substRate = substRate1.getValue();
        for (int [] i : ratemap1) {
        	int i1 = i[1];
        	int i2 = i[2];
            //relativeRates[i[0]] += rates1[i[1]][i[2]] * substRate;
            relativeRates[i[0]] += relrates1[i1*3+i2 + (i2 > i1 ? -1 : 0)] * substRate;
        }
	}


    @Override
    public double[] getFrequencies() {
    	if (frequencies != null) {
    		return frequencies.getFreqs();
    	}
    	
    	double [] freqs = new double[nrOfStates];
		double [] f1 = model1.getFrequencies();
		double [] f3 = model3.getFrequencies();
		double [] f2 = model2.getFrequencies();

		int i = 0;
		for (int i1 = 0; i1 < 4; i1++) {
			for (int j1 = 0; j1 < 4; j1++) {
				for (int k1 = 0; k1 < 4; k1++) {
					int k = index[i++];
					if (k >= 0) {
						freqs[k] += f1[i1] * f2[j1] * f3[k1];
					}
				}
			}
		}
    	return freqs;
    }

	@Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrix = true;
        return true;
    }

	@Override
    public boolean canHandleDataType(DataType dataType) {
    	return dataType instanceof Aminoacid;
    }

	public static void main(String[] args) {
		for (int i = 0; i < 750; i++) {
			System.out.println(i + " " + Math.exp(i));
		}
		
		
		ParameterisedAminoAcidModel m = new ParameterisedAminoAcidModel();
		m.initAndValidate();
		System.out.println(m.ratemap1.size());
		System.out.println(m.ratemap2.size());
		System.out.println(m.ratemap3.size());
	}
}
