package obama;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beastfx.app.tools.Application;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.util.LogFile;
import beastfx.app.util.XMLFile;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.inference.MCMC;
import beast.evolution.substitutionmodel.OBAMAModel;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.parser.XMLParser;

public class OBAMAAnalyser extends Runnable {
	public Input<XMLFile> xmlFileInput = new Input<>("xml", "XML file specifying the set of models", Validate.REQUIRED);
	public Input<LogFile> traceFileInput = new Input<>("log","trace log file containing output of a bModelTest analysis", Validate.REQUIRED);
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);

	NumberFormat formatter = new DecimalFormat("##0.00");     

	@Override
	public void initAndValidate() {
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
				    	modelName = modelName.replaceAll("OBAMA_", "");
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
				Log.info(modelCount[1] + "\t" + formatter.format(modelCount[1] * 100.0 / modelTrace.length) + "\tUse external frequencies");

			}
		}
		Log.info("Done");
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
		new Application(new OBAMAAnalyser(), "OBAMA Analyser", args);

	}

}
