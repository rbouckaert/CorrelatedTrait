package tsa.evolution.substitutionmodel;

import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.EigenSystem;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

public class MTJEigenSystem implements EigenSystem {

    private final int stateCount;

    public MTJEigenSystem(int stateCount) {

        this.stateCount = stateCount;

    }

    /**
     * set instantaneous rate matrix
     */
    public EigenDecomposition decomposeMatrix(double[][] AA) {

	        DenseMatrix DM = new DenseMatrix(AA);
	        EVD evd = new EVD(stateCount, true, true);
	        try {
				evd.factor(DM);
			} catch (NotConvergedException e2) {
				// TODO Auto-generated catch block
				e2.printStackTrace();
			}
	        DenseMatrix evec = evd.getLeftEigenvectors();
	        double [] Evec = ((DenseMatrix)evec.transpose()).getData();
	        double [] Eval = evd.getRealEigenvalues();
	        
//	        DenseMatrix I = Matrices.identity(n);
//	        DenseMatrix AI = I.copy();
//	        Matrix ivec = evec.solve(I, AI);
//	        Ievc = ((DenseMatrix)ivec).getData();
	        double [] Ievc = ((DenseMatrix)evd.getRightEigenvectors()).getData();
	        return new EigenDecomposition(Evec, Ievc, Eval);
	        
    }

}

