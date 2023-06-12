package tsa.correlatedcharacters.polycharacter;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;

@Description("A prior reporting on whether two traits are evolving dependently or independently"
		+ " in a CorrelatedSubstitutionModel.")
public class IndependencyLogger extends BEASTObject implements Loggable {
	public Input<CorrelatedSubstitutionModel> csmInput = new Input<CorrelatedSubstitutionModel>("model",
			"The CorrelatedSubstitutionModel this logger is reporting");
	
	protected Object trueOutput = true; 
	protected Object falseOutput = false; 
	
	@Override
	public void init(PrintStream out) {
		int components = csmInput.get().getShape().length;
		for (int component1 = 0; component1 < components; ++component1) {
			for (int component2 = 0; component2 < component1; ++component2) {
		        out.printf("%s_%d_depends_on_%d\t", getID(), component1, component2);				
		        out.printf("%s_%d_depends_on_%d\t", getID(), component2, component1);				
			}			
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		CorrelatedSubstitutionModel csm = csmInput.get();
		int components = csmInput.get().getShape().length;
		for (int component1 = 0; component1 < components; ++component1) {
			for (int component2 = 0; component2 < component1; ++component2) {
				if (csm.depends(component1, component2)) {
					out.print(1);
				} else {
					out.print(0);
				}
				out.print("\t");				
				if (csm.depends(component2, component1)) {
					out.print(1);
				} else {
					out.print(0);
				}
				out.print("\t");				
			}			
		}
		
	}

	@Override
	public void close(PrintStream out) {
        // nothing to do		
	}

	@Override
	public void initAndValidate() {}

}
