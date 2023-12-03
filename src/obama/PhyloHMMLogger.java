package obama;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.util.Randomizer;

@Description("Logs sample of HMM states of a Phylo-HMM")
public class PhyloHMMLogger extends BEASTObject implements Loggable {
	final public Input<PhyloHMM> phyloHMMInput = new Input<>("phyloHMM", "the phyloHMM for which to log its states", Validate.REQUIRED);
	
	
	PhyloHMM _phyloHMM;
	
	@Override
	public void initAndValidate() {
		_phyloHMM = phyloHMMInput.get();
	}

	@Override
	public void init(PrintStream out) {

		for (int i = 0; i < _phyloHMM.getSiteCount(); i++) {
			out.print((getID() != null ? getID() : "site") + i + "\t");  
		}
	}

	@Override
	public void log(long sample, PrintStream out) {

		_phyloHMM.calculateLogP();
		switch (_phyloHMM.hmmAlgorithmInput.get()) {
		case backwardforward:
			_phyloHMM.backward();
			break;

		case Viterbi:
			_phyloHMM.backwardViterbi();
			break;
		}
		
		for (int i = 0; i < _phyloHMM.getSiteCount(); i++) {
			out.print(Randomizer.randomChoicePDF(_phyloHMM.getHMMpartials()[i]) + "\t");  
		}
	}

	@Override
	public void close(PrintStream out) {

	}

}
