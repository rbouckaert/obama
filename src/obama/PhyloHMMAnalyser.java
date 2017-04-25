package obama;

import java.awt.Color;
import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.BEASTVersion2;
import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.app.util.Utils;
import beast.app.util.XMLFile;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.util.LogAnalyser;
import beast.util.XMLParser;
import beast.util.XMLParserException;

@Description("Visualises PhyloHMM log created by PhyloHMMLogger")
public class PhyloHMMAnalyser extends Runnable {
	public Input<LogFile> traceFileInput = new Input<>("file","trace log file containing output of a PhyloHMMLogger", Validate.REQUIRED);
	public Input<String> prefixInput = new Input<>("prefix", "prefix of the entry in the log file containing the sites traces (default 'site')" , "site");
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	public Input<Boolean> useBrowseInput = new Input<>("useBrowserForVisualisation", "use default web browser for visualising the output. "
			+ "Since not all browsers support all features, the alternative is to use an internal viewer, which requires an up to date Java 8 version.", true);
	public Input<XMLFile> xmlFileInput = new Input<>("xml", "XML file specifying the set of models. If specified, labels will be sourced from the XML file.");

	@Override
	public void initAndValidate() {
		int burnin = burninInput.get();
		if (burnin >= 100) {
			throw new IllegalArgumentException("burnin is a percentage and should not be larger than 100");
		}
	}

	LogAnalyser analyser;
	PhyloHMM phyloHMM = null;

	@Override
	public void run() throws Exception {
		File file = traceFileInput.get();
		String prefix = prefixInput.get();
		int burnin = burninInput.get();
		if (burnin < 0) {
			burnin = 0;
		}

		analyser = new LogAnalyser(file.getAbsolutePath(), burnin, false, false);
		
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
		svg.append("<svg width=\"" + (n + 100) + "\" height=\"" + (traceCount + 30 + stateCount * 25) + "\">\n");
		svg.append("<style  type=\"text/css\">\n");
		for (int i = 0; i < stateCount; i++) {
			Color color = new Color(Color.HSBtoRGB(((float) i)/stateCount, 0.9f, 0.9f));
			svg.append(".state" + i + " {fill:rgb(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ");"
					+ "stroke:rgb(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ")}\n");
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
		// labels
		String [] label = getLabels(stateCount);
		for (int i = 0; i < stateCount; i++) {
			int y = (45 + traceCount + 25 * i);
			svg.append("<rect x=\"10\" y=\"" + y + "\" width=\"50\" height=\"1\" class=\"state" + i + "\" "
					+ " style=\"stroke-width:3;\"/>\n");
			svg.append("<text x=\"100\" y=\"" + (y + 5) + "\" class=\"ticktext\" text-anchor=\"start\">" + label[i] +"</text>\n");
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

		String instance = "0";
		String jsPath = "/tmp";
		if (phyloHMM != null) {
			String dotty = toDotty(label);
			String dotty2 = dotty;
			dotty2 = dotty2.replaceAll("\n", "\\\\\n");
			dotty2 = dotty2.replaceAll("'", "&quot;");
			
			jsPath = getJavaScriptPath();
			String html = "<html>\n" + 
					"<title>BEAST " + new BEASTVersion2().getVersionString() + ": Phylo-HMM-Analyser</title>\n" +
					"<header>\n" + 
					"<link rel='stylesheet' type='text/css' href='css/style.css'>\n" +
					"<script data-main='" + jsPath + "/main' src='" + jsPath + "/requirejs/require.js'></script>\n" + 
					" \n" + 
					"<script>\n" + 
					"requirejs.config({\n" + 
					"    //By default load any module IDs from js-directory\n" + 
					"    baseUrl: '"+jsPath+"',\n" + 
					"    //except, if the module ID starts with 'app',\n" + 
					"    //load it from the js/app directory. paths\n" + 
					"    //config is relative to the baseUrl, and\n" + 
					"    //never includes a '.js' extension since\n" + 
					"    //the paths config could be for a directory.\n" + 
					"    paths: {\n" + 
					"        d3: 'd3/d3',\n" + 
					"        'dot-checker': 'graphviz-d3-renderer/dist/dot-checker',\n" + 
					"        'layout-worker': 'graphviz-d3-renderer/dist/layout-worker',\n" + 
					"        worker: 'requirejs-web-workers/src/worker',\n" + 
					"        renderer: 'graphviz-d3-renderer/dist/renderer'\n" + 
					"    }\n" + 
					"});\n" +
					"</script>\n" +
					"\n" +
					"</header>\n" +
					"<body>\n" +
					"\n" +
					"<svg id='graph' width='624' height='624'></svg>\n" +
					"<div id=\"img\" onclick=\"continueExecution()\">Create downloadable image</div>\n" + 
					"\n" + 
					"<script>\n" + 
					"  continueExecution = function() {\n" + 
					"var html = d3.select(\"svg\")\n" + 
					"        .attr(\"title\", \"bModelTest\")\n" + 
					"        .attr(\"version\", 1.1)\n" + 
					"        .attr(\"xmlns\", \"http://www.w3.org/2000/svg\")\n" + 
					"        .node().innerHTML;\n" + 
					"d3.select(\"#img\")\n" + 
					"        .html(\"Right-click on this preview and choose Save as\")\n" + 
					"        .append(\"img\")\n" + 
					"        .attr(\"src\", \"data:image/svg+xml;base64,\"+ btoa(html));\n" + 
					"}\n" + 
					"require(['renderer'],\n" + 
					"  function (renderer) {\n" + 
					"\n" + 
					//"var client = new XMLHttpRequest();\n" + 
					//"client.open('GET', 'bModelTest" + instance + ".dot');\n" + 
					//"client.onreadystatechange = function() {\n" + 
					//"  dotSource = client.responseText;\n" +
					" dotSource = '" + dotty2 + "'\n" +
					"  // initialize svg stage\n" + 
					"  renderer.init('#graph');\n" + 
					"\n" + 
					"  // update stage with new dot source\n" + 
					"  renderer.render(dotSource);\n" + 
					//"}\n" + 
					//"client.send();\n" + 
					"\n" + 
					"});\n" + 
					"</script>\n" + 
					"\n" + 
					svg.toString() +
					"</body>\n" +
					"</html>";
			try {
		        FileWriter outfile = new FileWriter(jsPath + "/PhyloHMM" + instance + ".html");
		        outfile.write(html);
		        outfile.close();
		        
		        outfile = new FileWriter("/tmp/p.dot");
		        outfile.write(dotty);
		        outfile.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
		
			jsPath = "/tmp";
			
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
		}
		
		if (useBrowseInput.get()) {
			try {
				openUrl("file://" + jsPath + "/PhyloHMM" + instance + ".html");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			new beast.app.tools.WebViewer("PhyloHMMAnalyser", "file://" + jsPath + "/PhyloHMM" + instance + ".html");
		}
	}

	private String getJavaScriptPath() {
		String classpath = System.getProperty("java.class.path");
		String[] classpathEntries = classpath.split(File.pathSeparator);
		String FILESEP = "/";
		if (Utils.isWindows()) {
			FILESEP = "\\\\";
		}
		for (String pathEntry : classpathEntries) {
			//Log.debug.print("Trying >" + pathEntry + "< ");
			if (new File(pathEntry).getName().toLowerCase().equals("bmodeltest.addon.jar")) {
				Log.debug.println("Got it!");
				File parentFile = (new File(pathEntry)).getParentFile().getParentFile();
				String parent = parentFile.getPath();
				return parent + FILESEP + "js";
			}
			//Log.debug.println("No luck ");
		}
		File f = new File(System.getProperty("user.dir") + FILESEP + ".." + FILESEP + "bModelTest" + FILESEP + "js");
		
		String jsPath = f.getPath();//System.getProperty("user.dir") + FILESEP + "js";
		jsPath = jsPath.replaceAll("\\\\", "/");
		//Log.debug.println("Using default: " + jsPath);
		return jsPath;
	}

	private String[] getLabels(int stateCount) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		if (xmlFileInput.get() != null && xmlFileInput.get().exists()) {
			XMLParser parser = new XMLParser();
			BEASTObject o = parser.parseFile(xmlFileInput.get());
			Set<Object> objects = new LinkedHashSet<>();
			getObjets(o, objects);
			for (Object o2 : objects) {
				if (o2 instanceof PhyloHMM) {
					phyloHMM = (PhyloHMM) o2;
					if (phyloHMM.stateLabelsInput.get() != null) {
						return phyloHMM.stateLabels;
					}
				}
			}
		}
		
		String [] label = new String[stateCount];
		for (int i = 0; i < stateCount; i++) {
			label[i] = "state " + i;
		}		
		return label;
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

	
	private String toDotty(String [] labels) {
		StringBuilder b = new StringBuilder();
		NumberFormat formatter = new DecimalFormat("##0.000");     

		b.append("digraph {\n");
			b.append(" graph [mindist=0.0, nodesep=0.25, ranksep=0.4]\n;");
			b.append(" node [fontsize=\"9\", style=\"solid\", color=\"#0000FF60\"];\n");
			b.append(" \n");
			for (int i = 0; i < labels.length; i++) {
				b.append(i + " [width=0.5, height=0.5, fixedsize=\"true\" color=\"#0000FF60\" label=\"" + labels[i] + "\"];\n");
			}
			b.append(" \n");
			
			Function rates = phyloHMM.rates;
			for (int i = 0; i < rates.getDimension(); i++) {
				double rate = getRate(rates, i);
				if (rate > 0) {
					int [] link = phyloHMM.getLink(i);
					int to = link[1];
					int from = link[0];
					b.append(from + " -> " + to + "[label=\"" + formatter.format(rate) + "\"];\n");
				}
			}
			
			b.append("\n");
			b.append("}\n");
		return b.toString();
	}

	private double getRate(Function rates, int i) {
		String id = ((BEASTInterface) rates).getID();
		List<String> labels = analyser.getLabels();
		for (int j = 0; j < labels.size(); j++) {
			if (labels.get(j).equals(id+i)) {
				Double [] trace = analyser.getTrace(id+i);
				double sum = 0;
				for (double d : trace) {
					sum += d;
				}
				double rate = sum / trace.length;
				return rate;
			}
		}
		// it is not in the trace file, so get it from the XML
		return rates.getArrayValue(i);
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
