package obama;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.app.util.XMLFile;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;
import beast.core.Input.Validate;
import beast.core.MCMC;
import beast.evolution.substitutionmodel.OBAMAModel;
import beast.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.util.LogAnalyser;
import beast.util.XMLParser;

public class BAMAAanalyser extends Runnable {
	public Input<XMLFile> xmlFileInput = new Input<>("xml", "XML file specifying the set of models", Validate.REQUIRED);
	public Input<LogFile> traceFileInput = new Input<>("log","trace log file containing output of a bModelTest analysis", Validate.REQUIRED);
	public Input<String> prefixInput = new Input<>("prefix", "prefix of the entry in the log file containing the substitution model trace (default 'modelIndex')" , "modelIndex");
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);

	NumberFormat formatter = new DecimalFormat("##0.00");     

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	@Override
	public void run() throws Exception {
		LogAnalyser trace = new LogAnalyser(traceFileInput.get().getPath(), burninInput.get(), false, false);
		
		XMLParser parser = new XMLParser();
		MCMC mcmc  = (MCMC) parser.parseFile(xmlFileInput.get());
		
		
		Set<Object> objects = new LinkedHashSet<>();
		getObjets(mcmc, objects);
		for (Object o : objects) {
			if (o instanceof OBAMAModel) {
				OBAMAModel bamaModel = (OBAMAModel) o;
				Log.info("Processsing " + bamaModel.getID() + "\n");

				List<EmpiricalSubstitutionModel> models = bamaModel.substModelInput.get();
				int [] modelCount = new int[models.size()];
				
				String logEntry = bamaModel.modelIndicatorInput.get().getID();
				if (!trace.getLabels().contains(logEntry) && logEntry.indexOf('.') > 0) {
					logEntry = logEntry.substring(0, logEntry.indexOf('.'));
				}
				if (trace.getLabels().contains(logEntry)) {
					Double [] modelTrace = trace.getTrace(logEntry);
					for (Double d : modelTrace) {
						modelCount[(int)(d+0.5)]++;
					}
					Log.info("count\tpercent\tmodel");
					Log.info("=============================");
					for (int i = 0; i < models.size(); i++) {
				    	EmpiricalSubstitutionModel model = models.get(i);
				    	String modelName = model.getID();
				    	if (modelName == null) {
				    		modelName = model.getClass().getSimpleName();
				    	}
				    	if (modelName != null && modelName.indexOf('.') > 0) {
				    		modelName = modelName.substring(0, modelName.lastIndexOf('.'));
				    	}
				    	modelName = modelName.replaceAll("BAMA_", "");
						Log.info(modelCount[i] + "\t" + formatter.format(modelCount[i] * 100.0 / modelTrace.length) + "\t" + modelName);
					}
				}

				// calc proportion of useExternalFreqs
				logEntry = bamaModel.useExternalFreqsInput.get().getID();
				if (!trace.getLabels().contains(logEntry) && logEntry.indexOf('.') > 0) {
					logEntry = logEntry.substring(0, logEntry.indexOf('.'));
				}
				modelCount = new int[2];
				Double [] modelTrace = trace.getTrace(logEntry);
				for (Double d : modelTrace) {
					modelCount[(int)(d+0.5)]++;
				}
				Log.info("\ncount\tpercent\tmodel");
				Log.info("=============================");
				Log.info(modelCount[0] + "\t" + formatter.format(modelCount[0] * 100.0 / modelTrace.length) + "\tUse model frequencies");
				Log.info(modelCount[1] + "\t" + formatter.format(modelCount[1] * 100.0 / modelTrace.length) + "\tUse empirical frequencies based on alignment");

			}
		}
	}

	private void getObjets(BEASTInterface bi, Set<Object> objects) {
		if (objects.contains(bi)) {
			return;
		}
		for (BEASTInterface bi2 : bi.listActiveBEASTObjects()) {
			getObjets(bi2, objects);
			objects.add(bi2);
		}
		
	}

	public static void main(String[] args) throws Exception {
		new Application(new BAMAAanalyser(), "BAMA Analyser", args);

	}

}
