package test.beast.evolution.operators;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.junit.Test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.RatesFreqsOperator;
import beast.evolution.substitutionmodel.EpochSubstitutionModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SVSGeneralSubstitutionModel;
import junit.framework.TestCase;

public class RatesFreqsOperatorTest extends TestCase {

	
	@Test
	public void testQ() {
		BooleanParameter mask = new BooleanParameter();
		mask.initByName("dimension", 20, "value", true);

		RealParameter rates = new RealParameter();
		rates.initByName("value", 
				       "1.0 1.0 1.0 1.0"
				+ " 1.0     1.0 1.0 1.0 "
				+ " 1.0 1.0     1.0 1.0 "
				+ " 1.0 1.0 1.0     1.0 "
				+ " 1.0 1.0 1.0 1.0 "
				+ "");
		
		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", "0.3 0.15 0.15 0.2 0.2");

		SVSGeneralSubstitutionModel Q = new SVSGeneralSubstitutionModel();
		Q.initByName(
//				"eigenSystem","beast.evolution.substitutionmodel.RobustEigenSystem",
				"frequencies", freqs, "rateIndicator", mask, "rates", rates, "symmetric", false);
		
		double [] matrix = new double[36];
		Q.getTransitionProbabilities(null, 100.0, 0.0, 1.0, matrix);
		printQ(matrix);
		
		double [] f = freqs.getFreqs();
		
		int n = 5;
		int [][] map = new int[n][n];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					map[i][j] = k++;
				}
			}			
		}
		int [] pair = new int[n*(n-1)];
		k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					pair[k++] = map[j][i];
				}
			}			
		}
		
		
		
		
		DecimalFormat formatter = new DecimalFormat("#.######");
		for (int i = 0; i < 20; i++) {
			double epsilon = 1e-1;
			double i0 = rates.getValue(i);
			double p0 = rates.getValue(pair[i]);
			rates.setValue(i, i0 + epsilon);
			rates.setValue(pair[i], p0 + epsilon);
			normalise(rates);
			double [] gradient = getGradient(matrix, Q, f);
			
			// should be all zeros
			System.out.print(i+" " + pair[i] +" [");
			for (k = 0; k < 5; k++) {
				System.out.print(formatter.format(gradient[k])+ ", ");
			}
			System.out.println("]");
		}

		
		RatesFreqsOperator operator = new RatesFreqsOperator();
		operator.initByName("rates", rates, "weight", 1.0);
		
		int N = 100000;
		for (int i = 0; i < N; i++) {
			operator.proposal();
		}
		double [] gradient = getGradient(matrix, Q, f);
		// should be all zeros
		System.out.print("diff in frequencies after " + N +" proposals = [");
		for (k = 0; k < 5; k++) {
			System.out.print(formatter.format(gradient[k])+ ", ");
			assertEquals(0.0, gradient[k], 1e-10);
		}
		System.out.println("]");
		
	}
	
	@Test
	public void testQ2() {
		BooleanParameter mask = new BooleanParameter();
		mask.initByName("dimension", 20, "value", true);

		RealParameter rates = new RealParameter();
		rates.initByName("value",
				  "         0.2408	0.0449	0.0379	0.0543	"
				+ "0.2403	        2.000	0.3149	0.0882	"
				+ "0.0708	2.000	        2.000	0.3263	"
				+ "0.0837	0.3813	2.0000	        1.6973	"
				+ "0.016	0.0438	0.1898	2.000"
				);

		
		rates.initByName("value",
			  "         0.2408	0.0449	0.0379	0.0543	"
			+ "0.2403	        2.2315	0.3149	0.0882	"
			+ "0.0708	2.2695	        4.5701	0.3263	"
			+ "0.0837	0.3813	4.8004	        1.6973	"
			+ "0.016	0.0438	0.1898	2.2983");
		
		rates.initByName("value", "1.0", "dimension", 20);
		
		rates.initByName("value",
				"     1  0  0  0" +
				"  1     2  0  0" +
				"  0  2     4  0" +
				"  0  0  4     8" +
				"  0  0  0  8   "
);
		jlabel = "UNI, rates 1 2 3 8, unequal frequencies";
		
////"				0.213		2.465E-4		1.139E-2	1.158E-2" + 	
////"	0.281					1.83			0.115		2.374E-2" + 	
////"	1.047E-4	5.227E-3					11.536		0.567	" + 
////"	0.236		4.179		3.746						0.822	" + 
////"	3.361E-3	7.807E-3	0.363			2.006				" + 			    
//				"");
		
		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", "0.453679653679654 0.296103896103896 0.14025974025974 0.071861471861472 0.038095238095238");
		
		SVSGeneralSubstitutionModel Q = new SVSGeneralSubstitutionModel();
		Q.initByName("frequencies", freqs, "rateIndicator", mask, "rates", rates, "symmetric", false);

		Frequencies freqs2 = new Frequencies();
		freqs2.initByName("frequencies", "0.753 0.163 0.084 0.00 0.0");
		SVSGeneralSubstitutionModel Q2 = new SVSGeneralSubstitutionModel();
		Q2.initByName("frequencies", freqs2, 
			"rateIndicator", 
			"     true true false false " +
			"true      true false false " +
			"true true      false false " +
			"false false false      false " +
			"false false false false     "
			, "rates", rates, "symmetric", false);

		EpochSubstitutionModel model = new EpochSubstitutionModel();
		model.initByName("epochDates", "6.0",
				"model", Q, "model", Q2, "frequencies", freqs);
		
		
		//double [] matrix = new double[36];
		double [][]matrix = new double[ITERATIONS][36];
		
		for (int i = 0; i < ITERATIONS; i++) {
//			model.getTransitionProbabilities(null, (i+1)*1.0, 0.0, 1.0, matrix);
//			System.out.println("model\nt = " + (i+1)*1);
//			printQ(matrix);

			Q.getTransitionProbabilities(null, (i)*DELTA, 0.0, 1.0, matrix[i]);
			System.out.println("Q\nt = " + (i)*DELTA);
			printQ(matrix[i]);

//			Q2.getTransitionProbabilities(null, (i+1)*1.0, 0.0, 1.0, matrix);
//			System.out.println("Q2\nt = " + (i+1)*1);
//			printQ(matrix);
		}
		display(matrix);
	}
	
	int ITERATIONS = 40;
	double DELTA = 0.25;
	String jlabel = "";


	private void display(double[][] matrix) {
		JPanel panel = new JPanel();
		panel.setLayout(new GridLayout(2, 3));
		for (int i = 0; i < 5; i++) {
			PDPanel graph = new PDPanel(i, matrix);
			panel.add(graph);
		}
		panel.add(new JLabel(jlabel));
		JFrame frame = new JFrame();
		frame.add(panel);
		frame.setSize(1024, 768);
		frame.setVisible(true);
		//JOptionPane.showMessageDialog(frame, "Done?"); 
		System.out.println("OK");
		
	}

	public static void main(String[] args) {
		RatesFreqsOperatorTest t = new RatesFreqsOperatorTest();
		t.testQ2();
	}
    /**
     * maps most significant digit to nr of ticks on graph *
     */
    final static int[] NR_OF_TICKS = new int[]{5, 10, 8, 6, 8, 10, 6, 7, 8, 9, 10};


    /* class for drawing information for a parametric distribution **/
    class PDPanel extends JPanel {
        // the length in pixels of a tick
        private static final int TICK_LENGTH = 5;

        // the right margin
        private static final int RIGHT_MARGIN = 20;

        // the margin to the left of y-labels
        private static final int MARGIN_LEFT_OF_Y_LABELS = 5;

        // the top margin
        private static final int TOP_MARGIN = 10;

        int m_nTicks;
        private static final long serialVersionUID = 1L;

        
        int src;
        double [][] matrix;
        
        PDPanel(int src, double [][] matrix) {
        	super();
        	this.src = src;
        	this.matrix = matrix;
        	
        }
        
        @Override
        public void paintComponent(java.awt.Graphics g) {

            Graphics2D g2d = (Graphics2D)g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            // Record current font, since drawError can take over part-way
            // through the call to drawGraph, which alters the graphics font size.
            Font originalFont = g.getFont();

            
            Color [] colors = new Color[] {Color.BLACK, Color.RED, Color.BLUE, Color.GREEN, Color.DARK_GRAY};
            	try {
            		for (int i = 0; i < 5; i++) {
                		drawGraph(50, g, i, colors[i]);
                	}
                } catch (Exception ex) {
                    System.out.println(ex.getMessage());
                    ex.printStackTrace();
                }
        }

        private void drawGraph(int labelOffset, Graphics g, int character, Color linecolor) {
            final int width = getWidth();
            final int height = getHeight();

            double minValue = 0;
            double maxValue = ITERATIONS*DELTA;
            double xRange = maxValue - minValue;
            // adjust yMax so that the ticks come out right
            double x0 = minValue;
            int k = 0;
            double f = xRange;
            double f2 = x0;
            while (f > 10) {
                f /= 10;
                f2 /= 10;
                k++;
            }
            while (f < 1 && f > 0) {
                f *= 10;
                f2 *= 10;
                k--;
            }
            f = Math.ceil(f);
            f2 = Math.floor(f2);
//			final int NR_OF_TICKS_X = NR_OF_TICKS[(int) f];
            for (int i = 0; i < k; i++) {
                f *= 10;
                f2 *= 10;
            }
            for (int i = k; i < 0; i++) {
                f /= 10;
                f2 /= 10;
            }
            //double adjXRange = f;

            xRange = xRange + minValue - f2;
            xRange = adjust(xRange);
            final int NR_OF_TICKS_X = m_nTicks;

            minValue = f2; //xRange = adjXRange;

            int points = ITERATIONS;
            int[] xPoints = new int[points];
            int[] yPoints = new int[points];
            double[] fyPoints = new double[points];
            double yMax = 1.0;
            for (int i = 0; i < points; i++) {
                //try {
                    fyPoints[i] = matrix[i][src * 5 + character]; 
                    		//getDensityForPlot(m_distr, minValue + (xRange * i) / points);
                //}
                if (Double.isInfinite(fyPoints[i]) || Double.isNaN(fyPoints[i])) {
                    fyPoints[i] = 0;
                }
                //fyPoints[i] = Math.exp(m_distr.logDensity(minValue + (xRange * i)/points));
                yMax = Math.max(yMax, fyPoints[i]);
            }

            yMax = 1.0;
            yMax = adjust(yMax);
            final int NR_OF_TICKS_Y = m_nTicks;

            // draw ticks on edge
            Font font = g.getFont();
            Font smallFont = new Font(font.getName(), font.getStyle(), font.getSize() * 2/3);
            g.setFont(smallFont);

            // collect the ylabels and the maximum label width in small font
            String[] ylabels = new String[NR_OF_TICKS_Y+1];
            int maxLabelWidth = 20;
            FontMetrics sfm = getFontMetrics(smallFont);
            for (int i = 0; i <= NR_OF_TICKS_Y; i++) {
                ylabels[i] = format(yMax * i / NR_OF_TICKS_Y);
                int stringWidth = sfm.stringWidth(ylabels[i]);
                if (stringWidth > maxLabelWidth) maxLabelWidth = stringWidth;
            }

            // collect the xlabels
            String[] xlabels = new String[NR_OF_TICKS_X+1];
            for (int i = 0; i <= NR_OF_TICKS_X; i++) {
                xlabels[i] = format(minValue + xRange * i / NR_OF_TICKS_X);
            }
            int maxLabelHeight = sfm.getMaxAscent()+sfm.getMaxDescent();

            maxLabelWidth = 20;
            maxLabelHeight = 10;
            
            int leftMargin = maxLabelWidth + TICK_LENGTH + 1 + MARGIN_LEFT_OF_Y_LABELS;
            int bottomMargin = maxLabelHeight + TICK_LENGTH + 1;

            int graphWidth = width - leftMargin - RIGHT_MARGIN;
            int graphHeight = height - TOP_MARGIN - bottomMargin - labelOffset;

            // DRAW GRAPH PAPER
            g.setColor(Color.WHITE);
//            g.fillRect(leftMargin, TOP_MARGIN, graphWidth, graphHeight);
            g.setColor(Color.BLACK);
            g.drawRect(leftMargin, TOP_MARGIN, graphWidth, graphHeight);

            for (int i = 0; i < points; i++) {
                xPoints[i] = leftMargin + graphWidth * i / points;
                yPoints[i] = 1 + (int) (TOP_MARGIN + graphHeight - graphHeight * fyPoints[i] / yMax);
            }
            g.setColor(linecolor);
            ((Graphics2D)g).setStroke(new BasicStroke(2f));
            g.drawPolyline(xPoints, yPoints, points);
            ((Graphics2D)g).setStroke(new BasicStroke(1f));

            for (int i = 0; i <= NR_OF_TICKS_X; i++) {
                int x = leftMargin + i * graphWidth / NR_OF_TICKS_X;
                g.drawLine(x, TOP_MARGIN + graphHeight, x, TOP_MARGIN + graphHeight + TICK_LENGTH);
                g.drawString(xlabels[i], x-sfm.stringWidth(xlabels[i])/2, TOP_MARGIN + graphHeight + TICK_LENGTH + 1 + sfm.getMaxAscent());
            }

            // draw the y labels and ticks
            for (int i = 0; i <= NR_OF_TICKS_Y; i++) {
                int y = TOP_MARGIN + graphHeight - i * graphHeight / NR_OF_TICKS_Y;
                g.drawLine(leftMargin - TICK_LENGTH, y, leftMargin, y);
                g.drawString(ylabels[i], leftMargin - TICK_LENGTH - 1 - sfm.stringWidth(ylabels[i]), y + 3);
            }

            int fontHeight = font.getSize() * 10 / 12;
            g.setFont(new Font(font.getName(), font.getStyle(), fontHeight));
        }
        
        private String format(double value) {
            StringWriter writer = new StringWriter();
            PrintWriter pw = new PrintWriter(writer);
            pw.printf("%.3g", value);
            if (value != 0.0 && Math.abs(value) / 1000 < 1e-320) { // 2e-6 = 2 * AbstractContinuousDistribution.solverAbsoluteAccuracy
            	pw.printf("*");
            }
            pw.flush();
            return writer.toString();
        }
        

        private double adjust(double yMax) {
            // adjust yMax so that the ticks come out right
            int k = 0;
            double y = yMax;
            while (y > 10) {
                y /= 10;
                k++;
            }
            while (y < 1 && y > 0) {
                y *= 10;
                k--;
            }
            y = Math.ceil(y);
            m_nTicks = NR_OF_TICKS[(int) y];
            for (int i = 0; i < k; i++) {
                y *= 10;
            }
            for (int i = k; i < 0; i++) {
                y /= 10;
            }
            return y;
        }
    }
	private void normalise(RealParameter rates) {
		double sum = 0;
		for (int i = 0; i < rates.getDimension(); i++) {
			sum += rates.getValue(i);
		}
		sum /= rates.getDimension();
		for (int i = 0; i < rates.getDimension(); i++) {
			rates.setValue(i, rates.getValue(i)/sum);
		}
		
		sum = 0;
		for (int i = 0; i < rates.getDimension(); i++) {
			sum += rates.getValue(i);
		}
		assertEquals(sum, rates.getDimension(), 1e-10);
	}

	private double [] getGradient(double [] matrix, SVSGeneralSubstitutionModel Q, double [] f) {
		Q.makeDirty();
		Q.getTransitionProbabilities(null, 100.0, 0.0, 1.0, matrix);
		double [] gradient = new double[5];
		for (int k = 0; k < 5; k++) {
			gradient[k] = f[k] - matrix[k]; 
		}
		Q.makeDirty();
		return gradient;
	}

	
	private void printQ(double [] matrix) {
		double [] d = new double[5];
		for (int i = 0; i < 5; i++) {
			System.arraycopy(matrix, i * 5, d, 0, 5);
			System.out.println(Arrays.toString(d));			
		}
	}
}
