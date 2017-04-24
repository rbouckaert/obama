package obama;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Param;

@Description("Specifies origin and target state for a rate in PhyloHMM")
public class Transition extends BEASTObject {
	private String from;
	private String to;
	
	public Transition(@Param(name="from", description="source state for origin/target pair") String from,
			@Param(name="to", description="target state for origin/target pair") String to) {
		this.from = from;
		this.to = to;
	}

	public String getFrom() {
		return from;
	}

	public void setFrom(String from) {
		this.from = from;
	}

	public String getTo() {
		return to;
	}

	public void setTo(String to) {
		this.to = to;
	}

	@Override
	public void initAndValidate() {}

}
