package obama;

import java.awt.Color;
import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Set;

import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.app.util.Utils;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.util.LogAnalyser;

@Description("Visualises PhyloHMM log created by PhyloHMMLogger")
public class PhyloHMMAnalyser extends Runnable {
	public Input<LogFile> traceFileInput = new Input<>("file","trace log file containing output of a bModelTest analysis", Validate.REQUIRED);
	public Input<String> prefixInput = new Input<>("prefix", "prefix of the entry in the log file containing the substitution model trace (default 'site')" , "site");
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	public Input<Boolean> useBrowseInput = new Input<>("useBrowserForVisualisation", "use default web browser for visualising the dot graph. "
			+ "Since not all browsers support all features, the alternative is to use an internal viewer, which requires an up to date Java 8 version.", true);

	@Override
	public void initAndValidate() {
		int burnin = burninInput.get();
		if (burnin >= 100) {
			throw new IllegalArgumentException("burnin is a percentage and should not be larger than 100");
		}
	}

	@Override
	public void run() throws Exception {
		File file = traceFileInput.get();
		String prefix = prefixInput.get();
		int burnin = burninInput.get();
		if (burnin < 0) {
			burnin = 0;
		}

		LogAnalyser analyser = new LogAnalyser(file.getAbsolutePath(), burnin, false, false);
		
		// sweep to determine nr of sites
		// and number of states
		double max = 0;
		int n = 0;
		for (String label : analyser.getLabels()) {
			if (label.startsWith(prefix)) {
				n++;
				for (Double d : analyser.getTrace(label)) {
					max = Math.max(d, max);
				}
			}
		}
		int stateCount = (int) (max + 1.5);
		int traceCount = analyser.getTrace(0).length;
		Log.warning(n + " sites with " + stateCount + " states " + traceCount + " samples");
		
		
		StringBuilder svg = new StringBuilder();
		svg.append("<svg width=\"" + (n + 100) + "\" height=\"" + (traceCount + 100) + "\">\n");
		svg.append("<style  type=\"text/css\">\n");
		for (int i = 0; i < stateCount; i++) {
			Color color = new Color(Color.HSBtoRGB(((float) i)/stateCount, 0.9f, 0.9f));
			svg.append(".state" + i + " {fill:rgb(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ")}\n");
		}
		svg.append(".ticktext {font: bold 8px}\n");
		svg.append("</style>\n");
		
		// bounding box
		svg.append("<rect x=\"9\" y=\"9\" width=\"" + (n + 2) + "\" height=\"" + (traceCount + 2) + "\" style=\"stroke-width:2;stroke:rgb(0,0,0)\"/>\n");

		// ticks
		for (int i = 10; i < n; i += 10) {
			if (i % 50 != 0) {
				svg.append("<line x1=\"" + (10 + i) + "\" y1=\"" + (10 + traceCount) + "\" x2=\"" + (10 + i) + "\" y2=\"" + (14 + traceCount) + "\" "
					+ " style=\"stroke-width:1;stroke:rgb(0,0,0)\"/>\n");
			}
		}
		for (int i = 50; i < n; i += 100) {
			svg.append("<line x1=\"" + (10 + i) + "\" y1=\"" + (10 + traceCount) + "\" x2=\"" + (10 + i) + "\" y2=\"" + (17 + traceCount) + "\" "
					+ " style=\"stroke-width:1;stroke:rgb(0,0,0)\"/>\n");
		}
		for (int i = 0; i < n; i += 100) {
			svg.append("<line x1=\"" + (10 + i) + "\" y1=\"" + (10 + traceCount) + "\" x2=\"" + (10 + i) + "\" y2=\"" + (20 + traceCount) + "\" "
					+ " style=\"stroke-width:1;stroke:rgb(0,0,0)\"/>\n");
			svg.append("<text x=\"" + (10 + i) + "\" y=\"" + (35 + traceCount) + "\" class=\"ticktext\" text-anchor=\"middle\">" + i +"</text>\n");
		}

		// content
		int x = 10;
		for (int i = 0; i < n; i++) {
			int [] count = new int[stateCount];
			for (Double d : analyser.getTrace(prefix + i)) {
				count[(int)(d + 0.5)]++;
			}
			int y = 10;
			for (int j = 0; j < stateCount; j++) {
				if (count[j] > 0) {
					svg.append("<rect x=\"" + x + "\" y=\"" + y + "\" width=\"1\" height=\"" + count[j] + "\" class=\"state" + j + "\"/>\n");
					y += count[j];
				}
			}
			x++;
		}
		svg.append("</svg>");

		
		String jsPath = "/tmp";
		String instance = "0";
		
		try {
	        FileWriter outfile = new FileWriter(jsPath + "/PhyloHMM" + instance + ".html");
	        outfile.write("<!DOCTYPE html>\n<html>\n<body>");
	        outfile.write(svg.toString());
	        outfile.write("</body>\n</html>");
	        outfile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if (useBrowseInput.get()) {
			try {
				openUrl("file://" + jsPath + "/PhyloHMM" + instance + ".html");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			new beast.app.tools.WebViewer("BModelAnalyser", "file://" + jsPath + "/PhyloHMM" + instance + ".html");
		}
	}

	void openUrl(String url) throws IOException {
		url = url.replaceAll(" ", "%20");
	    if(Desktop.isDesktopSupported()){
	        Desktop desktop = Desktop.getDesktop();
	        try {
	            desktop.browse(new URI(url));
	            return;
	        } catch (IOException | URISyntaxException e) {
	            // TODO Auto-generated catch block
	            e.printStackTrace();
	        }
	    }
	    if (Utils.isWindows()) {
	    	Runtime rt = Runtime.getRuntime();
	    	rt.exec( "rundll32 url.dll,FileProtocolHandler " + url);
	    } else if (Utils.isMac()) {
	    	Runtime rt = Runtime.getRuntime();
	    	rt.exec( "open" + url);
	    } else {
	    	// Linux:
	    	Runtime rt = Runtime.getRuntime();
	    	String[] browsers = {"epiphany", "firefox", "mozilla", "konqueror",
	    	                                 "netscape","opera","links","lynx"};

	    	StringBuffer cmd = new StringBuffer();
	    	for (int i=0; i<browsers.length; i++) {
	    	     cmd.append( (i==0  ? "" : " || " ) + browsers[i] +" \"" + url + "\" ");
	    	}
	    	rt.exec(new String[] { "sh", "-c", cmd.toString() });
	    }
	    
	   }
	
	public static void main(String[] args) throws Exception {
		new Application(new PhyloHMMAnalyser(), "Phylo-HMM Analyser", args);
	}
}
