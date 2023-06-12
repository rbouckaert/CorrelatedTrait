package tsa.correlatedcharacters.polycharacter;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import babel.tools.MatrixVisualiserBase;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beastfx.app.tools.LogAnalyser;

@Description("Visualises correlated subsitutition model rates matrix")
public class RateMatrixVisualiser extends MatrixVisualiserBase {
	final public Input<File> inFile = new Input<>("in", "trace file containing substitution matrix rates", new File("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("/tmp/matrix.svg"));

	@Override
	public void initAndValidate() {
	}

//	@Override
//	public void run() throws Exception {
// all the work is done in the base class
//	}

	@Override
	public double [][] getMatrix() {
		double [] r = new double[8];
		try {
			LogAnalyser tracelog = new LogAnalyser(inFile.get().getPath(), burnInPercentageInput.get(), true, false);
			CorrelatedSubstitutionModel s = new CorrelatedSubstitutionModel();
			
			final StringBuilder buf = new StringBuilder();
			PrintStream out = new PrintStream(new ByteArrayOutputStream(), restoreFromFile) {
				public void print(String s) {
					buf.append(s);
				};
			};
			s.init(out);
			String [] strs = buf.toString().split("\t");
			for (int i = 0; i < 8; i++) {
				r[i] = mean(tracelog.getTrace(strs[i]));
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.err.println(Arrays.toString(r));
		return new double[][] {
		  //"00","01","11", "10"
//			{0.0, r[0], 0.0, r[1]},
//		    {r[2], 0.0, r[3], 0.0},
//		    {0.0, r[6], 0.0, r[7]},
//		    {r[4], 0.0, r[5], 0.0},
		  // "11", "10", "00","01"
			{0.0, r[7], 0.0, r[6]},
		    {r[5], 0.0, r[4], 0.0},
		    {0.0, r[1], 0.0, r[0]},
		    {r[3], 0.0, r[2], 0.0},
		};
	}
	
	private double mean(Double[] trace) {
		double sum = 0;
		for (double d : trace) {
			sum += d;
		}
		return sum / trace.length;
	}

	@Override
	public String[] getLabels(double[][] rates) {
		return new String[] {"11", "10", "00","01"};
		// return new String[] {"00","01","11", "10"};
	}

	
	@Override
	public String getFileName() {
		return outputInput.get().getPath();
	}
	
	public static void main(String[] args) throws Exception {
		new Application(new RateMatrixVisualiser(), "Correlated rate matrix visualiser", args);
	}
}
