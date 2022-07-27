package obama;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.BooleanParameter;

@Description("Variable state size HMM variant of PhyloHMM")
public class PhyloHMMVar extends PhyloHMMSparse {
	final public Input<BooleanParameter> isActiveInput = new Input<>("isActive", "array of flags indicating which partition is active", Validate.REQUIRED);

	int maxHMMStateCount;
	
	BooleanParameter isActive;
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		isActive = isActiveInput.get();
		int sum = 0;
		for (int i = 0; i < isActive.getDimension(); i++) {
			sum += isActive.getArrayValue(i);
		}
		if (sum > HMMStateCount) {
			throw new IllegalArgumentException("number of active states in isActive (" + sum +") should be no more "
					+ "than the HMM state count (" + HMMStateCount +")"); 
		}
		maxHMMStateCount = HMMStateCount;
	}
	

	
	@Override
	protected boolean requiresRecalculation() {
		if (isActive.somethingIsDirty()) {
			int sum = 0;
			for (int i = 0; i < isActive.getDimension(); i++) {
				sum += isActive.getArrayValue(i);
			}
			if (sum > maxHMMStateCount) {
				resize(sum);
			}
		}
		return super.requiresRecalculation();
	}



	private void resize(int sum) {
		// TODO Auto-generated method stub
		
	}
}
