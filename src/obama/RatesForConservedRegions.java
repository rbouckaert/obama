package obama;

import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.parameter.RealParameter;

@Description("HMM rate matrix for Example 2 in Siepel, A. and Haussler, D., 2005. Phylogenetic hidden Markov models. In Statistical methods in molecular evolution (pp. 325-351). Springer New York.")
public class RatesForConservedRegions extends CalculationNode implements Function, Loggable {
	public Input<RealParameter> lambdaInput = new Input<>("lambda","auto correlation parameter lambda", Validate.REQUIRED);
	public Input<Integer> states = new Input<>("stateCount", "number of states k. Rates between states are (1.0-lambda)/k and rates to stay in states are lambda + (1.0-lambda)/k", 3);

	double [] rates;
	int k;
	RealParameter lambda;
	boolean needsUpdate = true;
	
	@Override
	public void initAndValidate() {
		k = states.get();
		rates = new double[k * k];
		lambda = lambdaInput.get();
		update();
	}

	private void update() {
		double b = (1.0 - lambda.getValue()) / k;
//		double a = lambda.getValue() + (k-1) * (1.0 - lambda.getValue()) / k;
		double a = lambda.getValue() + (1.0 - lambda.getValue()) / k;

		int u = 0;
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				if (i == j) {
					rates[u++] = a;
				} else {
					rates[u++] = b;
				}
			}
		}
		needsUpdate = false;
	}

	@Override
	public int getDimension() {
		return rates.length;
	}

	@Override
	public double getArrayValue() {
		if (needsUpdate) {
			update();
		}
		return rates[0];
	}

	@Override
	public double getArrayValue(int dim) {
		if (needsUpdate) {
			update();
		}
		return rates[dim];
	}
	
	@Override
	protected void store() {
		needsUpdate = true;
		super.store();
	}
	
	@Override
	protected void restore() {
		needsUpdate = true;
		super.restore();
	}

	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return super.requiresRecalculation();
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < rates.length; i++) {
			out.append(getID() + i + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < rates.length; i++) {
			out.append(rates[i] + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {		
	}
}
