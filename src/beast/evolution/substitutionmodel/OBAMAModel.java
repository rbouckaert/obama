package beast.evolution.substitutionmodel;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Node;

@Description("Substitution model that can average over a number of amino acid substitution models " +
		"as well as switch between the model's frequencies and external frequencies (as for example " +
		"empirical frequencies informed by an alignment).")
public class OBAMAModel extends GeneralSubstitutionModel {
	final public Input<BooleanParameter> useExternalFreqsInput = new Input<>("useExternalFreqs", "if false, use", new BooleanParameter("false"));
	final public Input<List<EmpiricalSubstitutionModel>> substModelInput = new Input<>("model", "empicial amino acid substitution model", new ArrayList<>(), Validate.REQUIRED);
	final public Input<IntegerParameter> modelIndicatorInput = new Input<>("modelIndicator", "index of the model in list of models that is used for its rates and frequencies", Validate.REQUIRED);

	BooleanParameter useExternalFreqs;
	IntegerParameter modelIndicator;
	List<EmpiricalSubstitutionModel> models;
	
	public OBAMAModel() {
		ratesInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() {
        frequencies = frequenciesInput.get();

		useExternalFreqs = useExternalFreqsInput.get();
		models = substModelInput.get();
		modelIndicator = modelIndicatorInput.get();
		if (modelIndicator.getUpper() > models.size() - 1) {
			Log.warning("Setting upper limit of " + modelIndicator.getID() + " to " + (models.size()-1) +".");
			modelIndicator.setUpper(models.size() - 1);
		}
		if (modelIndicator.getLower() < 0) {
			Log.warning("Setting lower limit of " + modelIndicator.getID() + " to 0.");
			modelIndicator.setLower(0);
		}
		
		
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;

        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e.getMessage());
		}
        //eigenSystem = new DefaultEigenSystem(m_nStates);

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates*(nrOfStates-1)];
        storedRelativeRates = new double[nrOfStates*(nrOfStates-1)];

	}
	
	
	@Override
    protected void setupRelativeRates() {
    	EmpiricalSubstitutionModel model = models.get(modelIndicator.getValue());
        System.arraycopy(model.m_empiricalRates, 0, relativeRates, 0, model.m_empiricalRates.length);
    }

	@Override
	public double[] getFrequencies() {
		if (useExternalFreqs.getValue()) {
			return super.getFrequencies();
		}
    	EmpiricalSubstitutionModel model = models.get(modelIndicator.getValue());
        return model.getFrequencies();
	}
	
	
    @Override
    public double[] getRateMatrix(Node node) {
    	EmpiricalSubstitutionModel model = models.get(modelIndicator.getValue());
        double[][] matrix = model.getEmpiricalRates();
        int states = matrix.length;
        double[] rates = new double[states * states];
        for (int i = 0; i < states; i++) {
            for (int j = i + 1; j < states; j++) {
                rates[i * states + j] = matrix[i][j];
                rates[j * states + i] = matrix[i][j];
            }
        }
        // determine diagonal
        for (int i = 0; i < states; i++) {
            double sum = 0;
            for (int j = i + 1; j < states; j++) {
                sum += rates[i * states + j];
            }
            rates[i * states + i] = -sum;
        }
        return rates;
    }	
	
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Aminoacid;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (useExternalFreqs.isDirtyCalculation()) {
			updateMatrix = true;
			return true;
		}
		if (modelIndicator.isDirtyCalculation()) {
			updateMatrix = true;
			return true;
		}
		return super.requiresRecalculation();
	}

}
