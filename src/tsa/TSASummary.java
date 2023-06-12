package tsa;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import beastfx.app.tools.Application;
import beastfx.app.util.LogFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Input.Validate;
import beastfx.app.tools.LogAnalyser;

@Description("Produces a summary for an analysis using the Correlated-character substitution model")
public class TSASummary extends Runnable {
	public Input<LogFile> traceFileInput = new Input<>("file","trace log file containing output of a TreeSetAnalyser analysis", Validate.REQUIRED);
	public Input<String> prefixInput = new Input<>("prefix", "prefix of the entry in the log file containing the substitution model trace (default 'rate_')" , "rate_");
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	public Input<String> iconInput = new Input<>("icons", "comma separated list to specify icons for first dimension, second dimension and both" , "");

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		File file = traceFileInput.get();
		String prefix = prefixInput.get();
		int burnin = burninInput.get();
		if (burnin < 0) {
			burnin = 0;
		}
		
		String [] icons = iconInput.get().split(",");
		
		LogAnalyser analyser = new LogAnalyser(file.getAbsolutePath(), burnin, false, true);
		
		double [] rates = new double[8];
		
		for (String label : analyser.getLabels()) {
			if (label.startsWith(prefix)) {
				String label0 = label.substring(prefix.length());
				switch (label0) {
				case "0,0->1,0": rates[0] = analyser.getMean(label); break;
				case "1,0->0,0": rates[1] = analyser.getMean(label); break;
				case "0,1->1,1": rates[2] = analyser.getMean(label); break;
				case "1,1->0,1": rates[3] = analyser.getMean(label); break;
				case "0,0->0,1": rates[4] = analyser.getMean(label); break;
				case "0,1->0,0": rates[5] = analyser.getMean(label); break;
				case "1,0->1,1": rates[6] = analyser.getMean(label); break;
				case "1,1->1,0": rates[7] = analyser.getMean(label); break;

				case "00->10": rates[0] = analyser.getMean(label); break;
				case "10->00": rates[1] = analyser.getMean(label); break;
				case "01->11": rates[2] = analyser.getMean(label); break;
				case "11->01": rates[3] = analyser.getMean(label); break;
				case "00->01": rates[4] = analyser.getMean(label); break;
				case "01->00": rates[5] = analyser.getMean(label); break;
				case "10->11": rates[6] = analyser.getMean(label); break;
				case "11->10": rates[7] = analyser.getMean(label); break;
				}
			}
		}
		
		String svg = 				 
				"<svg width='450' height='450'  xmlns:xlink='http://www.w3.org/1999/xlink' xmlns='http://www.w3.org/2000/svg'>\n" +
				 " <defs>\n" +
				 "    <radialGradient id='grad1' cx='50%' cy='50%' r='50%' fx='50%' fy='50%'>\n" +
				 "      <stop offset='0%' style='stop-color:rgb(255,255,255);\n" +
				 "      stop-opacity:0'></stop>\n" +
				 "      <stop offset='100%' style='stop-color:rgb(0,0,255);stop-opacity:1'></stop>\n" +
				 "    </radialGradient>\n" +
				 "\n" +
				 "    <marker id='head' viewBox='0 0 20 20' refX='2' refY='10' markerUnits='strokeWidth' markerWidth='8' markerHeight='6' orient='auto'>\n" +
				 "        <path d='M 0 0 L 20 10 L 0 20 z' fill='black'/>\n" +
				 "      </marker>\n" +
				 "    </defs>    \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[4])+"' fill='none' stroke='black' d='M80,25 Q200,0 290,25'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[5])+"' fill='none' stroke='black' d='M318,75 Q200,100 110,75'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[6])+"' fill='none' stroke='black' d='M80,325 Q200,300 290,325'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[7])+"' fill='none' stroke='black' d='M318,375 Q200,400 110,375'></path>  \n" +
				 "\n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[0])+"' fill='none' stroke='black' d='M25,80 Q0,200 25,290'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[1])+"' fill='none' stroke='black' d='M75,318 Q100,200 75,110'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[2])+"' fill='none' stroke='black' d='M325,80 Q300,200 325,290'></path>  \n" +
				 "  <path marker-end='url(#head)' stroke-width='" + arrowWidth(rates[3])+"' fill='none' stroke='black' d='M375,318 Q400,200 375,110'></path>  \n" +
				 "\n" + 
			     (icons.length < 3 ?
				 "   <circle cx='50' cy='50' r='40' stroke='black' stroke-width='4' fill='url(#grad1)' />\n" + 
				 "   <circle cx='350' cy='50' r='40' stroke='black' stroke-width='4' fill='url(#grad1)' />\n" + 
				 "   <circle cx='50' cy='350' r='40' stroke='black' stroke-width='4' fill='url(#grad1)' />\n" + 
				 "   <circle cx='350' cy='350' r='40' stroke='black' stroke-width='4' fill='url(#grad1)' />\n"
				 : "") +	
				 "\n" + 
				 "	<text x='25' y='62' font-family='Verdana' font-size='35'>00</text>\n" + 
				 "	<text x='325' y='62' font-family='Verdana' font-size='35'>01</text>\n" + 
				 "	<text x='25' y='362' font-family='Verdana' font-size='35'>10</text>\n" + 
				 "	<text x='325' y='362' font-family='Verdana' font-size='35'>11</text>\n" + 
				 "\n" + 
				 "	<text transform='translate(25,200)rotate(90)' text-anchor='middle' font-size='35'>" + format(rates[0])+ "</text>\n" + 
				 "	<text transform='translate(100,200)rotate(90)' text-anchor='middle' font-size='35'>" + format(rates[1])+ "</text>\n" + 
				 "	<text transform='translate(275,200)rotate(90)' text-anchor='middle' font-size='35'>" + format(rates[2])+ "</text>\n" + 
				 "	<text transform='translate(350,200)rotate(90)' text-anchor='middle' font-size='35'>" + format(rates[3])+ "</text>\n" + 
				 "	<text x='200' y='50' text-anchor='middle' font-size='35'>" + format(rates[4])+ "</text>\n" + 
				 "	<text x='200' y='125' text-anchor='middle' font-size='35'>" + format(rates[5])+ "</text>\n" + 
				 "	<text x='200' y='300' text-anchor='middle' font-size='35'>" + format(rates[6])+ "</text>\n" + 
				 "	<text x='200' y='375' text-anchor='middle' font-size='35'>" + format(rates[7])+ "</text>\n" + 
			     (icons.length == 3 ?
					"  <image x='0' y='302' width='100' height='100' xlink:href='" + icons[0] + "'/>\n" +
					"  <image x='305' y='2' width='100' height='100' xlink:href='" + icons[1] + "'/>\n" +
					"  <image x='305' y='302' width='100' height='100' xlink:href='" + icons[2] + "'/>\n"
					: ""					
				 ) +
				 "\n" + 
				 "Sorry, your browser does not support inline SVG.\n" + 
				 "</svg> \n";

		

		
		try {
			File tmpFile0 = File.createTempFile("TSASummary", ".svg");
			FileWriter outfile = new FileWriter(tmpFile0);
	        outfile.write(svg);
	        outfile.close();

			String html = "<!DOCTYPE html>\n" + 
					 "<html>\n" + 
					 "<body>\n" + 
					 "\n" + 
					 "<h1>TreeSetAnalyser "+ file.getName() +"</h1>\n" + 
					 "\n" +
					 svg +
					 " \n" +
					 "<p><a href='file://" + tmpFile0.getPath() + "'>svg file</a>" +
					 "</body>\n" + 
					 "</html>\n";

			File tmpFile = File.createTempFile("TSASummary", ".html");
	        outfile = new FileWriter(tmpFile);
	        outfile.write(html);
	        outfile.close();
			
			Application.openUrl("file://" + tmpFile.getPath());

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	double arrowWidth(double rate) {
		double width = 30 * rate / 8.0;
		if (width < 0.5) {
			width = 0.5;
		}
		return width;
	}
	
	String format(double d) {
		NumberFormat formatter = new DecimalFormat("##0.00");     
		return formatter.format(d);
	}

	public static void main(String[] args) throws Exception {
		new Application(new TSASummary(), "TSA Summary", args);
	}

}
