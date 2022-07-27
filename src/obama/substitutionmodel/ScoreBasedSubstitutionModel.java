package obama.substitutionmodel;


import java.lang.reflect.InvocationTargetException;

import beast.base.core.Description;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

@Description("Amino acid substituion model based on scores.")
abstract public class ScoreBasedSubstitutionModel extends GeneralSubstitutionModel {
	
	public ScoreBasedSubstitutionModel() {
		ratesInput.setRule(Validate.OPTIONAL);
	}

	double [] Q;
	
	@Override
	public void initAndValidate() {		
        frequencies = frequenciesInput.get();

        updateMatrix = true;
        nrOfStates = 20;

        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e.getMessage());
		}

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[20 * 19];
        storedRelativeRates = new double[relativeRates.length];

        Q = new double[20*19];
		setUpQMatrix();
	}
	
	private void setUpQMatrix() {
		int [][]scores = getScores();
		double [] freqs = getFrequencies();
		
		// normalise to log probs
		// PMatrix[i][j] = exp(score[i][j])/(f[i] f[j])
		double [][] PMatrix = new double[20][20];
		for (int i = 0; i < 20; i++) {
			for (int j = i; j < 20; j++) {
				PMatrix[i][j] = Math.exp(scores[i][j])/(freqs[i] * freqs[j]);
				PMatrix[i][j] = PMatrix[i][j];
			}			
		}
		
		// make sure probabilities add to 1
		for (int i = 0; i < 20; i++) {
			double sum = 0;
			for (int j = 0; j < 20; j++) {
				sum += PMatrix[i][j];
			}			
			for (int j = 0; j < 20; j++) {
				PMatrix[i][j] /= sum;
			}			
		}
		
		eigenDecomposition = eigenSystem.decomposeMatrix(PMatrix);
		
		// take the log of the matrix
        int i, j, k;
        double temp;
        double[] iexp = new double[nrOfStates * nrOfStates];
        // Eigen vectors
        double[] Evec = eigenDecomposition.getEigenVectors();
        // inverse Eigen vectors
        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
        // Eigen values
        double[] Eval = eigenDecomposition.getEigenValues();
        for (i = 0; i < nrOfStates; i++) {
            temp = Math.log(Eval[i]);
            for (j = 0; j < nrOfStates; j++) {
                iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
            }
        }

        int u = 0;
        for (i = 0; i < nrOfStates; i++) {
            for (j = 0; j < nrOfStates; j++) {
            	if (i != j) {
	                temp = 0.0;
	                for (k = 0; k < nrOfStates; k++) {
	                    temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
	                }
	
	                Q[u] = Math.abs(temp);
	                u++;
            	}
            }
        }
	}

	@Override
	public void setupRelativeRates() {
    	setUpQMatrix();
    	System.arraycopy(Q, 0, relativeRates, 0, Q.length);
    }

    /** return scores matrix **/
    abstract int [][] getScores();
	

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Aminoacid;
	}

}
