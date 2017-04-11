import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class Magnetic {
	
	private final double mu=1., dx = 1.;//, dt = 0.1;
	private int N=100;
	private int[] iPlus = new int[N+1], iMinus = new int[N+1];
	private double[][] A_z = new double[N+1][N+1];
	private double[][] j_z = new double[N+1][N+1];
	private double[] B_ijk = new double[3];
	
	public Magnetic(int N, double A_z0) {
		this.N = N;
		this.init();
	}
	
	private void setPlusMinus(){
		iPlus = new int[N+1];
		iMinus = new int[N+1];
		for(int i=1; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
	}
	
	private void init(){
		
		this.setPlusMinus();
		A_z = new double[N+1][N+1];
		j_z = new double[N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++){
					if (i==N/2 && j==N/2) j_z[i][j] = 1.0;
					else j_z[i][j] = 0.;
				}
		
		this.setNew(1.);
		
	}
	
	private void getB_ijk(int i, int j){
		B_ijk[0] = (A_z[i][iPlus[j]]-A_z[i][iMinus[j]])/(2.*dx); 
		B_ijk[1] = -(A_z[iPlus[i]][j]-A_z[iMinus[i]][j])/(2.*dx);
		B_ijk[2] = 0.;
	}
	
	private double setNew(double omega){
		
		double A_z_ij_n = 0.;
		double error = 0.;
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++){
				A_z_ij_n = A_z[i][j];
				A_z[i][j] = (A_z[iPlus[i]][j]+
							A_z[iMinus[i]][j]+A_z[i][iPlus[j]]+
							A_z[i][iMinus[j]]+
							j_z[i][j]*mu)/4.0;
				A_z[i][j] = (1.-omega)*A_z_ij_n + omega * A_z[i][j];
				error += Math.abs(A_z_ij_n - A_z[i][j]);
			}
		return error;
	}
	
	public int evolve(String A_zFile, double omega, boolean write) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		double error = 0.;
		int count = 0;
		while(true){
			
			error = this.setNew(omega);
			count++;
			if(count%100==0 & count > 1000){
				if(error < 1e-4)
					break;
			}
		}
		
		if(write){
			PrintWriter writer1 = new PrintWriter(A_zFile, "UTF-8");
			for(int i=1; i<N; i++){
				for(int j=1; j<N; j++){
					this.getB_ijk(i, j);
					writer1.println(i + " " + j +  " " + 
								j_z[i][j] + " " + A_z[i][j] + " " +
								B_ijk[0] + " " + B_ijk[1] + " " + B_ijk[2]);
				}
				writer1.println();
			}
			writer1.close();
		}
		return count;
		
	}
	
	public void SOR(String SORfile) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		PrintWriter writer1 = new PrintWriter(SORfile, "UTF-8");
		int Nomega = 50;
		double omega, d_omega = 1./(double)(Nomega);
		int nt = 0;
		for(int w=0; w<Nomega+5; w++){
			omega = 1. + d_omega * w;
			nt = this.evolve("", omega, false);
			writer1.println(omega + " " + nt);
			this.init();
		}
		writer1.close();
	}
	
	public static void main(String[] args) 
			throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		Magnetic wire = new Magnetic(100, 0.0);
		//wire.evolve("MagneticWithCutoff.dat", 1., true);
		wire.SOR("magneticSOR_1e-4.dat");
	}
	
}
